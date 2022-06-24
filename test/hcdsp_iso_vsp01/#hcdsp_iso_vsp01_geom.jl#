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

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz.seis"));

# get geometry
sx = SeisMain.ExtractHeader(hzx,"sx");
sy = SeisMain.ExtractHeader(hzx,"sy");
gl = SeisMain.ExtractHeader(hzx,"gelev");

dt = ext_zz.d1;
nt = ext_zz.n1;

ntr = ext_zz.n2;

drz = 16.6 ;       # m 
rz_init = 1350.0f0;  # m
rz_end  = 1850.01f0; # m
rz_axis = range(rz_init, rz_end, step=drz);
nr = length(rz_axis)

nt = 218;     # number of time samples (int64)
ns = 205;     # number of sources within lines
nsline = 205; # number of source lines

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;

@assert ns*nsline*nr == ntr

ozx = zeros(eltype(dzx),nt,ns,nsline,nr);
ozy = zeros(eltype(dzy),nt,ns,nsline,nr);
ozz = zeros(eltype(dzz),nt,ns,nsline,nr);
#count = zeros(Int,ns,nsline,nr);

for i in eachindex(hzx)

    # i-th trace header
    h = hzx[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    ozx[:,sxi,syi,gel] += dzx[:,i]
    ozy[:,sxi,syi,gel] += dzy[:,i]
    ozz[:,sxi,syi,gel] += dzz[:,i]

    #    count[sxi,syi,gel] += 1;
end

# select a crg
icrg = 15;
ozx = ozx[:,:,:,icrg];
ozy = ozy[:,:,:,icrg];
ozz = ozz[:,:,:,icrg];

#nshots
nshots = ns*nsline

dt = 0.012f0

# record length
rec_length = dt*nt

# total time for conventional acqusition
tconv = rec_length*nshots

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
#β = 0.9;

# firing times
#tau = Float32.(sort(rand(nshots).* β .* tconv));

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

# Inverse crime: observed data
b  = bFwd(ozz);

# Pseudo-deblended
db = bAdj(b);

# Patching
psize = nextpow.(2,(100,20,20));
polap = (10,20,20);
 smin = ntuple(x->1,Val(ndims(db)));
 smax = (nt,ns,nsline);

dpatch,pid = fwdPatchOp(ozz,psize,polap,smin,smax);

db2 = adjPatchOp(dpatch,pid,psize,polap,smin,smax);

# Threshold schedule
@everywhere Pi, Pf, K = 99.9, 0.1, 101
sched = HCDSP.thresh_sched(dpatch,K,Pi,Pf,"exp") ./ 5;

# Define a patching operator
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

# Step-size
α = Float32(0.1);

# Deblending by inversion
pgd_fkt,it_pgd_fkt = pgdls!(bFwd, bAdj, b, zero(db),
                            proj!, (psize,polap,smin,smax,sched);
                            ideal = ozz, α=α,
                            verbose=true,
                            maxIter=K,
                            ε=Float32(1e-16));
# Reg param
λ = Float32(0.5);

# Deblending by inversion    
fp_fkt,it_fp_fkt = red_fp!(bFwd, bAdj, b, zero(db), λ,
                           proj!, (psize,polap,smin,smax,sched);
                           ideal = ozz,
                           verbose=true,
                           max_iter_o=K,
                           max_iter_i=10,
                           ε=Float32(1e-16));

# red reg param
λ = Float32(0.5);
γ = Float32(0.5);

# Deblending by inversion    
admm_fkt,it_admm_fkt = red_admm!(bFwd, bAdj, b, zero(db), λ, γ,
                                 proj!, (psize,polap,smin,smax,sched);
                                 ideal = ozz,
                                 verbose=true,
                                 max_iter_o=K,
                                 max_iter_i1=10,
                                 max_iter_i2=1,
                                 ε=Float32(1e-16));

j=50;
tmp = [ozz[:,:,j] db[:,:,j] pgd_fkt[:,:,j] fp_fkt[:,:,j] admm_fkt[:,:,j]];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=90,vmin=-a,vmax=a,cmap="gray")

#=
# Define a patching operator
function ssa_proj!(state, (psize, polap, smin, smax, dt, fmin, fmax, rank))
    # output allocation
    out = copy(state.x)

    # get iteration:
    it = state.it;

    # define ssa function
    fssa(δ) = fx_process(δ,dt,fmin,fmax,fast_ssa_lanc,(rank))

    # apply patching on input
    patches,pid = fwdPatchOp(out, psize, polap, smin, smax);
    
    # fk_thresh all patches
    patches .= pmap(fssa,patches);
    
    # Rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);

    # Return
    return out
end

fmin=0;fmax=45;rank=3;

# Deblending by inversion
pgd_ssa,it_pgd_ssa = pgdls!(bFwd, bAdj, b, zero(db),
                            ssa_proj!, (psize,polap,smin,smax,dt,fmin,fmax,rank);
                            ideal = ozz, α=α,
                            verbose=true,
                            maxIter=K,
                            ε=Float32(1e-6));
=#
