pwd()

using Distributed
addprocs(19);

@everywhere dev_dir = "/dev/Breno_GOM/projects";
cd(dev_dir)
pwd()

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()

@everywhere using Revise
@everywhere using LinearAlgebra
@everywhere using FFTW
@everywhere using HCDSP

using PyPlot
using SeisMain, SeisPlot

# data dir home
dir_path  = joinpath(dev_dir,"files/iso_vsp01")
dx_path = joinpath(dir_path,"iso_vsp01_zx.seis");
dy_path = joinpath(dir_path,"iso_vsp01_zy.seis");
dz_path = joinpath(dir_path,"iso_vsp01_zz.seis");

# read data
dzx,hzx,ext_zx = SeisRead(dx_path);
dzy,hzy,ext_zy = SeisRead(dy_path);
dzz,hzz,ext_zz = SeisRead(dz_path);

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

ozx = zeros(eltype(dzx),nt,ns,nsline,nr);
ozy = zeros(eltype(dzy),nt,ns,nsline,nr);
ozz = zeros(eltype(dzz),nt,ns,nsline,nr);

for i in eachindex(hzz)

    # i-th trace header
    h = hzz[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    ozx[:,sxi,syi,gel] += dzx[:,i]
    ozy[:,sxi,syi,gel] += dzy[:,i]
    ozz[:,sxi,syi,gel] += dzz[:,i]
end

#nshots
nshots = ns*nsline

# sampling time
dt = 0.012

# record length
rec_length = dt*nt

# shooting grid in indexes
grid = [(ix,iy) for ix in 1:ns, iy in 1:nsline]

# boat shooting positions
boat1 = copy(grid[:,1:103]) # from (1,1)
boat2 = reverse(reverse(grid[:,104:end],dims=1),dims=2) # from (205,205)

# number of shots per boat
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
tau = Vector{Float64}(undef,nshots); tau .= tmp;

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

bFwd(x) = QBlendOp(x, PARAM, "fwd");
bAdj(x) = QBlendOp(x, PARAM, "adj");

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
    tmp = PadOp(out; nin=n, npad=npad, flag="fwd")
    
    # fft
    tmp .= fft(tmp)

    # threshold
    threshold!(tmp,sched)
    
    # Truncate
    out .= PadOp( ifft(tmp); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

x_pgd_fkt = size(ozx) |> zeros;
y_pgd_fkt = size(ozy) |> zeros;
z_pgd_fkt = size(ozz) |> zeros;

x_fp_fkt = size(ozx) |> zeros;
y_fp_fkt = size(ozy) |> zeros;
z_fp_fkt = size(ozz) |> zeros;

x_admm_fkt = size(ozx) |> zeros;
y_admm_fkt = size(ozy) |> zeros;
z_admm_fkt = size(ozz) |> zeros;

it_pgd_fkt_snr  = Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef,nr);
it_fp_fkt_snr   = Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float64}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float64}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float64}}(undef,nr);

# tolerance
ε=1e-6;

for icrg in 1:nr

    println("Deblend $icrg-th CRG")
    
    # Select a CRG
    x = Float64.(ozx[:,:,:,icrg]);
    y = Float64.(ozy[:,:,:,icrg]);
    z = Float64.(ozz[:,:,:,icrg]);
    q = quaternion(x,y,z);

    # Inverse crime: model observed data
    b  = bFwd(q);

    # Pseudo-deblended
    db = bAdj(b);

    # Patching
    dpatch,_ = fwdPatchOp(db,psize,polap,smin,smax);

    # Threshold scheduler (it should be √.(sched) but that is too much)
    sched = HCDSP.thresh_sched(dpatch,K,Pi,Pf,"exp") ./ 5;
   
    # RED reg param
    λ = 0.5;

    # Deblending by inversion    
    tmp,tmp_it = red_fp!(bFwd, bAdj, b, zero(db), λ,
                         proj!, (psize,polap,smin,smax,sched);
                         ideal = q,
                         verbose=true,
                         max_iter_o=K,
                         max_iter_i=10,
                         ε=ε);
    
    # Store inversion results
    x_fp_fkt[:,:,:,icrg] .= imagi.(tmp);
    y_fp_fkt[:,:,:,icrg] .= imagj.(tmp);
    z_fp_fkt[:,:,:,icrg] .= imagk.(tmp);

    it_fp_fkt_snr[icrg] = tmp_it[:snr];
    it_fp_fkt_mis[icrg] = tmp_it[:misfit];
 
    # RED reg & admm param
    λ = 0.5; γ = 0.5;

    # Deblending by inversion    
    tmp,tmp_it = red_admm!(bFwd, bAdj, b, zero(db), λ, γ,
                           proj!, (psize,polap,smin,smax,sched);
                           ideal = q,
                           verbose=true,
                           max_iter_o=K,
                           max_iter_i1=10,
                           max_iter_i2=1,
                           ε=ε);

    # Store inversion results
    x_admm_fkt[:,:,:,icrg] .= imagi.(tmp);
    y_admm_fkt[:,:,:,icrg] .= imagj.(tmp);
    z_admm_fkt[:,:,:,icrg] .= imagk.(tmp);

    it_admm_fkt_snr[icrg] = tmp_it[:snr];
    it_admm_fkt_mis[icrg] = tmp_it[:misfit];

    # PGD step-size (< 1/β ≈ 0.5)
    α = 0.1;

    # Deblending by inversion
    tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                        proj!, (psize,polap,smin,smax,sched);
                        ideal = q, α=α,
                        verbose=true,
                        maxIter=K,
                        ε=ε);
    
    # Store inversion results
    x_pgd_fkt[:,:,:,icrg] .= imagi.(tmp);
    y_pgd_fkt[:,:,:,icrg] .= imagj.(tmp);
    z_pgd_fkt[:,:,:,icrg] .= imagk.(tmp);

    it_pgd_fkt_snr[icrg] = tmp_it[:snr];
    it_pgd_fkt_mis[icrg] = tmp_it[:misfit];

end

SeisWrite("files/iso_vsp01/iso_vsp01_zz_pgd_qfkt.seis", pgd_fkt,hzz,ext)
SeisWrite("files/iso_vsp01/iso_vsp01_zz_fp_qfkt.seis",  fp_fkt,hzz,ext)
SeisWrite("files/iso_vsp01/iso_vsp01_zz_admm_qfkt.seis",admm_fkt,hzz,ext)

fname = joinpath(dir_path,"all_iter_hist_qfkt.h5")

using HDF5

fid = h5open(fname,"w")

create_group(fid,"admm/snr")
create_group(fid,"fp/snr")
create_group(fid,"pgd/snr")

create_group(fid,"admm/mis")
create_group(fid,"fp/mis")
create_group(fid,"pgd/mis")

for i in 1:nr
    fid["admm"]["snr"]["$i"] = it_admm_fkt_snr[i]
    fid["fp"]["snr"]["$i"]   = it_fp_fkt_snr[i]
    fid["pgd"]["snr"]["$i"]  = it_pgd_fkt_snr[i]

    fid["admm"]["mis"]["$i"] = it_admm_fkt_mis[i]
    fid["fp"]["mis"]["$i"]   = it_fp_fkt_mis[i]
    fid["pgd"]["mis"]["$i"]  = it_pgd_fkt_mis[i]
end

close(fid)
