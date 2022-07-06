pwd()

using Distributed
addprocs(9);

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

@everywhere include(joinpath(dev_dir,"dev/deblending/src/deblend_module.jl"))
@everywhere import Main.deblend_module: BlendOp

using PyPlot
using SeisMain, SeisPlot

# data dir home
dir_path  = joinpath(dev_dir,"files/iso_vsp01")
data_path = joinpath(dir_path,"iso_vsp01_zz.seis");

# read data
dzz,hzz,ext_zz = SeisRead(data_path);

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

ozz = zeros(eltype(dzz),nt,ns,nsline,nr);
pgd_fkt  = similar(ozz);
fp_fkt   = similar(ozz);
admm_fkt = similar(ozz);

for i in eachindex(hzz)

    # i-th trace header
    h = hzz[i]

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

# parameter selection
@everywhere fmin=0;fmax=50;rank=5;

# Params for patching
psize = nextpow.(2,(100,20,20));
polap = (10,20,20);
smin = (1,1,1);
smax = (nt,ns,nsline);

# Define a patching operator
function proj!(state, (psize, polap, smin, smax, dt, fmin, fmax, rank))
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

it_pgd_fkt_snr  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_snr   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Float32}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float32}}(undef,nr);

# tolerance
ε=Float32(1e-4);

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

    # PGD step-size (< 1/β ≈ 0.5)
    α = Float32(0.2);

    # Deblending by inversion
    tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                        proj!, (psize,polap,smin,smax,dt,fmin,fmax,rank);
                        ideal = tmpz, α=α,
                        verbose=true,
                        maxIter=K,
                        ε=ε);  
    
    # Store inversion results
    pgd_fkt[:,:,:,icrg] .= tmp;
    it_pgd_fkt_snr[icrg] = tmp_it[:snr];
    it_pgd_fkt_mis[icrg] = tmp_it[:misfit];
    
    # RED reg param
    λ = Float32(0.2);

    # Deblending by inversion    
    tmp,tmp_it = red_fp!(bFwd, bAdj, b, zero(db), λ,
                         proj!, (psize,polap,smin,smax,dt,fmin,fmax,rank);
                         ideal = tmpz,
                         verbose=true,
                         max_iter_o=K,
                         max_iter_i=10,
                         ε=ε);
    
    # Store inversion results
    fp_fkt[:,:,:,icrg] .= tmp;
    it_fp_fkt_snr[icrg] = tmp_it[:snr];
    it_fp_fkt_mis[icrg] = tmp_it[:misfit];
 
    # RED reg & admm param
    λ = Float32(0.1);
    γ = Float32(0.1);

    # Deblending by inversion    
    tmp,tmp_it = red_admm!(bFwd, bAdj, b, zero(db), λ, γ,
                           proj!, (psize,polap,smin,smax,dt,fmin,fmax,rank);
                           ideal = tmpz,
                           verbose=true,
                           max_iter_o=K,
                           max_iter_i1=10,
                           max_iter_i2=1,
                           ε=ε);

    # Store inversion results
    admm_fkt[:,:,:,icrg] .= tmp;
    it_admm_fkt_snr[icrg] = tmp_it[:snr];
    it_admm_fkt_mis[icrg] = tmp_it[:misfit];
end

ext = ext_zz;
ext.n2 = 205; ext.n3=205; ext.n4=31;
ext.label2="source line"; ext.unit2="index";
ext.label3="source number"; ext.unit3="index";
ext.label4="receiver"; ext.unit4="index";

SeisWrite("files/iso_vsp01/iso_vsp01_zz_pgd_ssa.seis",pgd_fkt,hzz,ext)
SeisWrite("files/iso_vsp01/iso_vsp01_zz_fp_ssa.seis",fp_fkt,hzz,ext)
SeisWrite("files/iso_vsp01/iso_vsp01_zz_admm_ssa.seis",admm_fkt,hzz,ext)

fname = joinpath(dir_path,"all_iter_hist_ssa.h5")

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
