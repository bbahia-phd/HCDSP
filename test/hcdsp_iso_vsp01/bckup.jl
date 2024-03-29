cd(joinpath(homedir(),"projects"))
pwd()

using Distributed

addprocs(Sys.CPU_THREADS-1);

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(homedir(),"projects/HCDSP"))

@everywhere using Revise
@everywhere using LinearAlgebra
@everywhere using FFTW
@everywhere using HCDSP

@everywhere include("./dev/deblending/src/deblend_module.jl")
@everywhere using Main.deblend_module

using PyPlot
using SeisMain, SeisPlot

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";

# read data
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz.seis"));

# get geometry
sx = SeisMain.ExtractHeader(hzz,"sx");
sy = SeisMain.ExtractHeader(hzz,"sy");
gl = SeisMain.ExtractHeader(hzz,"gelev");

nt = 218;     # number of time samples (int64)
ns = 205;     # number of sources within lines
nsline = 205; # number of source lines
nr = 31;      # number of receivers down-well
drz = 16.6 ;  # source and receiver spacing (m )

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;

@assert ns*nsline*nr == ntr

ozz = zeros(eltype(dzz),nt,ns,nsline,nr);
pgd_fkt  = similar(ozz);
fp_fkt   = similar(ozz);
admm_fkt = similar(ozz);

for i in eachindex(hzx)

    # i-th trace header
    h = hzx[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    ozz[:,sxi,syi,gel] += dzz[:,i]

end

#nshots
nshots = ns*nsline

dt = 0.012f0

# record length
rec_length = dt*nt

grid = [(ix,iy) for ix in 1:ns, iy in 1:nsline]

# boat shooting positions
boat1 = copy(grid[:,1:103])
boat2 = reverse(reverse(grid[:,104:end],dims=1),dims=2)

nb1 = prod(size(boat1))
nb2 = prod(size(boat2))

# boat shooting times
tb1 = collect(0 : rec_length : rec_length*(nb1-1))

tmp = round.(Int, (tb1 .+ 3 .* rand(nb1)) ./ dt)
tb2 = sort(tmp .* dt)[1:nb2]

i=10
[tb1[i] tb2[i]]

nsx = [];
nsy = [];
tmp = [];

for i in eachindex(boat1)
    if i <= nb2 
        append!(tmp,[tb1[i];tb2[i]])
        append!(nsx,[boat1[i][1];boat2[i][1]])
        append!(nsy,[boat1[i][2];boat2[i][2]])
    else
        append!(tmp,tb1[i])
        append!(nsx,boat1[i][1])
        append!(nsy,boat1[i][2])
    end
end

isx = Vector{Int64}(undef,nshots); isx .= nsx;
isy = Vector{Int64}(undef,nshots); isy .= nsy;
tau = Vector{Float32}(undef,nshots); tau .= tmp;

# blending factor
t_conv = rec_length*nshots
t_blend = maximum(tau) + rec_length -1
β = t_conv/t_blend;

#
PARAM = (nt = nt,
         nx = ns,
         ny = nsline,
         dt = dt,
         tau = tau,
         sx = isx,
         sy = isy)

bFwd(x) = BlendOp(x, PARAM, "fwd");
bAdj(x) = BlendOp(x, PARAM, "adj");

# Threshold schedule
@everywhere Pi, Pf, K = 99.9, 0.1, 101

# Params for patching
psize = nextpow.(2,(100,20,20));
polap = (10,20,20);
smin = (1,1,1);
smax = (nt,ns,nsline);

# Define a patched projection operator
function proj!(state, (psize, polap, smin, smax, sched))
    # output allocation
    out = copy(state.x)

    # get iteration:
    it = state.it;

    # apply patching on input
    patches,pid = fwdPatchOp(out, psize, polap, smin, smax);
    
    # fk_thresh all patches
    patches .= pmap(fk_thresh,
                    patches,
                    repeat([sched[it]], length(patches)));
    
    # Rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);

    # Return
    return out
end

@everywhere function fk_thresh(IN::AbstractArray,sched::AbstractFloat)

    out = copy(IN)
    
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad
    tmp = complex.(PadOp(out; nin=n, npad=npad, flag="fwd"))
    
    # fft
    fft!(tmp)

    # threshold
    threshold!(tmp,sched)
    
    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

it_pgd_fkt_snr  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_snr   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Float32}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float32}}(undef,nr);

for icrg in 1:nr

    println("Deblend $icrg-th CRG")
    
    # Select a CRG
    tmpz = ozz[:,:,:,icrg];

    # Inverse crime: model observed data
    b  = bFwd(tmpz);

    # Pseudo-deblended
    db = bAdj(b);

    # Patching
    dpatch,_ = fwdPatchOp(db,psize,polap,smin,smax);

    # Threshold scheduler
    sched = HCDSP.thresh_sched(dpatch,K,Pi,Pf,"exp") ./ 5;

    # PGD step-size (< 1/β ≈ 0.5)
    α = Float32(0.1);

    # Deblending by inversion
    tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                        proj!, (psize,polap,smin,smax,sched);
                        ideal = tmpz, α=α,
                        verbose=true,
                        maxIter=K,
                        ε=Float32(1e-16));
    
    # Store inversion results
    pgd_fkt[:,:,:,icrg] .= tmp;
    it_pgd_fkt_snr[icrg] = tmp_it[:snr];
    it_pgd_fkt_mis[icrg] = tmp_it[:misfit];
    
    # RED reg param
    λ = Float32(0.5);

    # Deblending by inversion    
    tmp,tmp_it = red_fp!(bFwd, bAdj, b, zero(db), λ,
                         proj!, (psize,polap,smin,smax,sched);
                         ideal = tmpz,
                         verbose=true,
                         max_iter_o=K,
                         max_iter_i=10,
                         ε=Float32(1e-16));
    
    # Store inversion results
    fp_fkt[:,:,:,icrg] .= tmp;
    it_fp_fkt_snr[icrg] = tmp_it[:snr];
    it_fp_fkt_mis[icrg] = tmp_it[:misfit];
 
    # RED reg & admm param
    λ = Float32(0.5);
    γ = Float32(0.5);

    # Deblending by inversion    
    tmp,tmp_it = red_admm!(bFwd, bAdj, b, zero(db), λ, γ,
                           proj!, (psize,polap,smin,smax,sched);
                           ideal = tmpz,
                           verbose=true,
                           max_iter_o=K,
                           max_iter_i1=10,
                           max_iter_i2=1,
                           ε=Float32(1e-16));

    # Store inversion results
    admm_fkt[:,:,:,icrg] .= tmp;
    it_admm_fkt_snr[icrg] = tmp_it[:snr];
    it_admm_fkt_mis[icrg] = tmp_it[:misfit];
end

j=50;
tmp = [ozz[:,:,j] db[:,:,j] pgd_fkt[:,:,j] fp_fkt[:,:,j]];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=90,vmin=-a,vmax=a,cmap="gray")
