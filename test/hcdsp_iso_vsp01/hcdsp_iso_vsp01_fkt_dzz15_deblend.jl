pwd()
cd(joinpath(homedir(),"projects"))

using Distributed

addprocs(4);

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(homedir(),"projects/HCDSP"))

@everywhere using Revise
@everywhere using LinearAlgebra
@everywhere using FFTW
@everywhere using HCDSP
@everywhere using Random

@everywhere include("./dev/deblending/src/deblend_module.jl")
@everywhere using Main.deblend_module

using PyPlot
using SeisMain, SeisPlot
using HDF5

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

ozz = zeros(eltype(dzz),nt,ns,nsline,nr);

for i in eachindex(hzz)

    # i-th trace header
    h = hzz[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    ozz[:,sxi,syi,gel] += dzz[:,i]

end

# Get one CRG
ozz = ozz[:,:,:,15];

#nshots
nshots = ns*nsline

dt = 0.012f0

# record length
rec_length = dt*nt

grid = [(ix,iy) for ix in 1:ns, iy in 1:nsline];

# boat shooting positions
boat1 = copy(grid[:,1:103]);
boat2 = reverse(reverse(grid[:,104:end],dims=1),dims=2);

nb1 = prod(size(boat1));
nb2 = prod(size(boat2));

# Random seed
Random.seed!(1992)

# boat shooting times
tb1 = collect(0 : rec_length : rec_length*(nb1-1));

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
@everywhere Pi, Pf, K = 99.9, 0.01, 101

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

nr = 1;
 
it_pgd_fkt_snr  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_snr   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Float32}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float32}}(undef,nr);

pgd_fkt  = similar(ozz);
fp_fkt   = similar(ozz);
admm_fkt = similar(ozz);

# Select a CRG
tmpz = ozz;

# Inverse crime: model observed data
b  = bFwd(tmpz);

# Pseudo-deblended
db = bAdj(b);

# Patching
dpatch,_ = fwdPatchOp(db,psize,polap,smin,smax);

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,K,Pi,Pf,"exp") ./ 5;

####################################

# PGD step-size (< 1/β ≈ 0.5)
α = Float32(0.5);

# Deblending by inversion
tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = tmpz, α=α,
                    verbose=true,
                    maxIter=K,
                    ε=Float32(1e-16));

# Store inversion results
pgd_fkt .= tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

####################################

# RED reg param
λ = Float32(0.2);

# Deblending by inversion    
tmp,tmp_it = red_fp!(bFwd, bAdj, b, zero(db), λ,
                     proj!, (psize,polap,smin,smax,sched);
                     ideal = tmpz,
                     verbose=true,
                     max_iter_o=K,
                     max_iter_i=10,
                     ε=Float32(1e-16));

# Store inversion results
fp_fkt .= tmp;
it_fp_fkt_snr = tmp_it[:snr];
it_fp_fkt_mis = tmp_it[:misfit];

# RED reg & admm param
λ = Float32(0.2); γ = Float32(0.2);

####################################

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
admm_fkt .= tmp;
it_admm_fkt_snr = tmp_it[:snr];
it_admm_fkt_mis = tmp_it[:misfit];

####################################

j=50;
tmp = [ozz[:,:,j] db[:,:,j] pgd_fkt[:,:,j] fp_fkt[:,:,j] admm_fkt[:,:,j]];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=90,vmin=-a,vmax=a,cmap="gray")

j=50;
tmp = [ozz[:,:,j] db[:,:,j] 5 .* (ozz .- pgd_fkt)[:,:,j] 5 .* (ozz .- fp_fkt)[:,:,j]  5 .* (ozz .- admm_fkt)[:,:,j]];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=99,vmin=-a,vmax=a,cmap="gray"

####################################
           
# data dir home
data_path = "/media/bbahia/DATA/seismic_data/iso_vsp01";

ext = ext_zz;
ext.n2 = 205; ext.n3 = 205; ext.n4 = 1;
ext.label2="source line";   ext.unit2="index";
ext.label3="source number"; ext.unit3="index";
ext.label4="receiver";      ext.unit4="index";

SeisWrite(joinpath(data_path,"iso_vsp01_zz_pgd_fkt.seis"), pgd_fkt,hzz,ext)
SeisWrite(joinpath(data_path,"iso_vsp01_zz_fp_fkt.seis"),  fp_fkt,hzz,ext)
SeisWrite(joinpath(data_path,"iso_vsp01_zz_admm_fkt.seis"),admm_fkt,hzz,ext)

####################################

fname = joinpath(data_path,"all_iter_hist_fkt.h5")

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

####################################

### read rec gains from qfkt
fname = joinpath(data_path,"qfkt/all_iter_hist_qfkt.h5")
fid = h5open(fname,"r")

q_admm_mis = fid["admm"]["mis"] |> read;
q_fp_mis = fid["fp"]["mis"]     |> read;
q_pgd_mis = fid["pgd"]["mis"]   |> read;

q_admm_snr = fid["admm"]["snr"] |> read;
q_fp_snr = fid["fp"]["snr"]     |> read;
q_pgd_snr = fid["pgd"]["snr"]   |> read;
close(fid)

####################################

q_fp_snr[q_fp_snr .== 0.0] .= NaN;

## read rec gains from qssa
fname = joinpath(data_path,"ssa/all_iter_hist_ssa.h5")
fid = h5open(fname,"r")

ssa_admm_mis = fid["admm"]["mis"]["15"] |> read;
ssa_fp_mis = fid["fp"]["mis"]["15"]     |> read;
ssa_pgd_mis = fid["pgd"]["mis"]["15"]   |> read;

ssa_admm_snr = fid["admm"]["snr"]["15"] |> read;
ssa_fp_snr = fid["fp"]["snr"]["15"]     |> read;
ssa_pgd_snr = fid["pgd"]["snr"]["15"]   |> read;
ssa_pgd_snr = vcat(ssa_pgd_snr,NaN .* ones(K-length(ssa_pgd_snr)))
close(fid)

####################################

## read rec gains from qssa
fname = joinpath(data_path,"qssa/all_iter_hist_qssa.h5")
fid = h5open(fname,"r")

qssa_admm_mis = fid["admm"]["mis"] |> read;
qssa_fp_mis = fid["fp"]["mis"]     |> read;
qssa_pgd_mis = fid["pgd"]["mis"]   |> read;

qssa_admm_snr = fid["admm"]["snr"] |> read;
qssa_fp_snr = fid["fp"]["snr"]     |> read;
qssa_pgd_snr = fid["pgd"]["snr"]   |> read;

qssa_pgd_snr[qssa_pgd_snr .== 0.0] .= NaN;

close(fid)

sfkt_snr = [it_pgd_fkt_snr it_fp_fkt_snr it_admm_fkt_snr];
qfkt_snr = [q_pgd_snr[15,:,3] q_fp_snr[15,:,3] q_admm_snr[15,:,3]];
sssa_snr = [ssa_pgd_snr ssa_fp_snr ssa_admm_snr];
qssa_snr = [qssa_pgd_snr[15,:,3] qssa_fp_snr[15,:,3] qssa_admm_snr[15,:,3]];

####################################

figname="scalar_vs_vector_plot";
close("all")
figure(figname,figsize=(20,4) .* 0.75)

subplot(141)
plot(1:K,sfkt_snr[:,1], label="PGD (FKT)", color="black")
plot(1:K,sfkt_snr[:,2], label="RED-FP (FKT)", color="red")
plot(1:K,sfkt_snr[:,3], label="RED-ADMM (FKT,"*L" I=1)", color="green")
title("(a)")
xlabel("Iteration number (k)",fontsize=15)
ylabel(L"R_z = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_2 }{\parallel \hat{\bf d}_k - {\bf d}^{o} \parallel^2_2} \right)}",fontsize=15)
ylim(0,100)
legend()

subplot(142)
plot(1:K,qfkt_snr[:,1], label="PGD (QFKT)", color="black")
plot(1:K,qfkt_snr[:,2], label="RED-FP (QFKT)", color="red")
plot(1:K,qfkt_snr[:,3], label="RED-ADMM (QFKT,"*L" I = 1)", color="green")
title("(b)")
xlabel("Iteration number (k)",fontsize=15)
ylabel(L"R_z = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_2 }{\parallel \hat{\bf d}_k - {\bf d}^{o} \parallel^2_2} \right)}",fontsize=15)
ylim(0,100)
legend()

subplot(143)
plot(1:K,sssa_snr[:,1], label="PGD (SSA)",color="black")
plot(1:K,sssa_snr[:,2], label="RED-FP (SSA)",color="red")
plot(1:K,sssa_snr[:,3], label="RED-ADMM (SSA,"*L" I = 1)",color="green")
title("(c)")
xlabel("Iteration number (k)",fontsize=15)
ylabel(L"R_z = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_2 }{\parallel \hat{\bf d}_k - {\bf d}^{o} \parallel^2_2} \right)}",fontsize=15)
ylim(0,60)
legend()

subplot(144)
plot(1:K,qssa_snr[:,1], label="PGD (QSSA)",color="black")
plot(1:K,qssa_snr[:,2], label="RED-FP (QSSA)",color="red")
plot(1:K,qssa_snr[:,3], label="RED-ADMM (QSSA,"*L" I = 1)",color="green")
title("(d)")
xlabel("Iteration number (k)",fontsize=15)
ylabel(L"R_z = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_2 }{\parallel \hat{\bf d}_k - {\bf d}^{o} \parallel^2_2} \right)}",fontsize=15)
ylim(0,60)
legend()

tight_layout()

savefig(joinpath("/home/bbahia/projects/files",figname))
