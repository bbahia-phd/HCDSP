
import Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using Revise

using LazyGrids
using HCDSP
using PyPlot
using SeisPlot, SeisMain
using StatsBase, Statistics
using LinearAlgebra
using FFTW

# kill me, kill me now, but this will give us access to wigb plotting 
using MATLAB
mat"addpath('./../../../QSEIS/src/Plotting/')"

function get_mode_irregular_data(x1,x2;ot=0.0, dt=0.004, nt=200, tau=[0.3],
                                p1=[0.0001],p2=[0.0],p3=[0.0],p4=[0.0],
                                amp=[1.0], f0=20.0)

    params_zx = (ot=ot, dt=dt, nt=nt, tau=tau,
    p1=p1, p2=p2, amp=amp, f0=f0)
    p = SeisLinearIrregularEvents3D(x1,x2; params_zx...);

    params_zy = (ot=ot, dt=dt, nt=nt, tau=tau,
    p1=p1, p2=p2, amp=amp, f0=f0)
    sv = SeisLinearIrregularEvents3D(x1,x2; params_zy...);

    params_zz = (ot=ot, dt=dt, nt=nt, tau=tau,
    p1=p1, p2=p2,amp=amp, f0=f0)
    sh = SeisLinearIrregularEvents3D(x1,x2; params_zz...);

    return (p,sv,sh)

end

function unmix(p,sv,sh)

    A = inv([0.75 0.15 0.1; 0.15 0.75 0.1; 0.1 0.15 0.75]);
    
    o1,o2,o3 = similar(p),similar(p),similar(p)
    for i in eachindex(p)
        tmp = A*[p[i]; sv[i]; sh[i]]
        o1[i] = tmp[1];
        o2[i] = tmp[2];
        o3[i] = tmp[3];
    end

    return o1,o2,o3
end

function mix(p,sv,sh)

    o1,o2,o3 = similar(p),similar(p),similar(p)
    for i in eachindex(p)
        o1[i] = 0.75*p[i] + 0.15*sv[i] + 0.1*sh[i];
        o2[i] = 0.15*p[i] + 0.75*sv[i] + 0.1*sh[i];
        o3[i] = 0.1*p[i] + 0.15*sv[i] + 0.75*sh[i];
    end

    return o1,o2,o3
end

# Ok. This is how Rongzhi is doing:
nx,ny,nt=30,30,200;
dx,dy,dt=20,20,0.004;

# midpoints
mx = collect(0:nx-1).*dx;
my = collect(0:ny-1).*dy;

(gmx,gmy) = ndgrid(mx,my);

jx = (rand(size(gmx)...) .* 2 .- 1) .* dx
jy = (rand(size(gmy)...) .* 2 .- 1) .* dy

jmx = gmx .+ jx; jmx[1,:] .= abs.(jmx[1,:]);
jmy = gmy .+ jy; jmy[:,1] .= abs.(jmy[:,1]);

px = [1/4000,  1/9000, -1/4000];
py = [1/4000, -1/4000, -1/6000];

τ  = [0.1,0.3,0.6];
amp = [1,-1,1];

f0 = 20.0;

# clean & pure seismic modes
p,sv,sh = get_mode_irregular_data(gmx,gmy);

# mixed observed displacements
dzz,dzy,dzx = mix(p,sv,sh);

# ideal
d0 = copy(quaternion(dzx,dzy,dzz));

# clean & pure irregular seismic modes
pi,svi,shi = get_mode_irregular_data(jmx,jmy);

# mixed observed displacements
diz,diy,dix = mix(pi,svi,shi);

# Plot using wigb from MATLAB
t = (0:nt-1)*dt
x = jmx[:,1]

# see irregular data (can plot anything like this tho)
mat"wigb($(diz[:,1,:]),1,$x,$t)"

# decimate
perc1=50;

# indexes
indx0=1:nx*ny                           # all indexes
indx1=HCDSP.decimate_indexes(dzz,perc1) # missing indexes
indx2=setdiff(indx0,indx1)              # observed indexes

# (random) decimating regular data (data regular observed = dro)
drz = reshape(dzz,nt,:);
drz[:,indx1] .= zero(eltype(drz));

dry = reshape(dzy,nt,:);
dry[:,indx1] .= zero(eltype(dry));

drx = reshape(dzx,nt,:);
drx[:,indx1] .= zero(eltype(drx));

# (random) decimating irregular data (data irregular observerd = dio)
diz = reshape(diz,nt,:);
diz[:,indx1] .= zero(eltype(diz));

diy = reshape(diy,nt,:);
diy[:,indx1] .= zero(eltype(diy));

dix = reshape(dix,nt,:);
dix[:,indx1] .= zero(eltype(dix));
 
# Sampling operator (no need)
T = SamplingOp(dix);

# extract the real irregular without binning
dix2 = dix[:,indx2];
diy2 = diy[:,indx2];
diz2 = diz[:,indx2];

# how is this binning dio2: dio2 is a bunch of (observed) traces together in matrix form (albeit 3D)
# You are placing dio2 into regular bins in the dio array.
# This would be the binning step, since dio was gridded accordingly.
dix[:,indx2] .= dix2;
diy[:,indx2] .= diy2;
diz[:,indx2] .= diz2;

# The 3D versions of the regular and irregular decimated data
dr3dx = reshape(drx,nt,nx,ny);
dr3dy = reshape(dry,nt,nx,ny);
dr3dz = reshape(drz,nt,nx,ny);

di3dx = reshape(dix,nt,nx,ny);
di3dy = reshape(diy,nt,nx,ny);
di3dz = reshape(diz,nt,nx,ny);

# 3D version of sampling
T = reshape(T,nt,nx,ny);

# vectorize the grid and get the coordinate (irregular) of the observations
jmx_obs = reshape(jmx,:,1)[indx2];
jmy_obs = reshape(jmy,:,1)[indx2];

# create the arrays containing the regular and irregular grids.
# These are input to the interpolator_kaiser_sinc_3d_freq
# Regular grid is 3D
h = zeros(nx,ny,2);
h[:,:,1] .= gmx;
h[:,:,2] .= gmy;

# Irregular grid is 2D. I think data will be too? Yes
htt = zeros(length(jmx_obs),2);
htt[:,1] = jmx_obs;
htt[:,2] = jmy_obs;

din = quaternion(dix2,diy2,diz2);

##########################################################################
# Define a patching operator
function proj!(state, (dt,fmin,fmax,rank))
    # output allocation
    out = copy(state.x)
    
    # get iteration:
    it = state.it;
    
    # fk_thresh
    out .= fx_process(out,dt,fmin,fmax,fast_qssa_lanc,(rank))
    
    # Return
    return out
end

fwd(x) = interp_ks3d(x,htt,h,3,10,"fwd")
adj(x) = interp_ks3d(x,htt,h,3,10,"adj")

dadj = adj(din);

# Step-size selection
α = 0.1;

# Number iterations
K = 101;

# f-x process
dt=0.004; fmin=0; fmax=80; rank=5;

# Reconstruction via PGD+SSA (I-MSSA)
out_ssa,it_ssa = pgdls!(fwd, adj, din, zero(dadj),
                        proj!, (dt,fmin,fmax,rank),
                        ideal=d0,
                        α = α, verbose=true,
                        maxIter=K, ε=1e-6);

L=length(it_ssa[:snr])
pgd_qssa_snrx = zeros(L);
pgd_qssa_snry = zeros(L);
pgd_qssa_snrz = zeros(L);

for i in 1:length(it_ssa[:snr])
    pgd_qssa_snrx[i]=it_ssa[:snr][i][1]
    pgd_qssa_snry[i]=it_ssa[:snr][i][2]
    pgd_qssa_snrz[i]=it_ssa[:snr][i][3]
end


# Reg param
λ = 8.0;

# Reconstruction via RED(FP)+SSA
red_ssa,red_it_ssa = red_fp!(fwd, adj, din, zero(dadj), λ,
                             proj!, (dt,fmin,fmax,rank),
                             ideal=d0,
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-6);


M=length(red_it_ssa[:snr])
fp_qssa_snrx = zeros(M);
fp_qssa_snry = zeros(M);
fp_qssa_snrz = zeros(M);

for i in 1:M
    fp_qssa_snrx[i]=red_it_ssa[:snr][i][1]
    fp_qssa_snry[i]=red_it_ssa[:snr][i][2]
    fp_qssa_snrz[i]=red_it_ssa[:snr][i][3]
end


# Reg param
λ = 8.0;

# ADMM-param
β = 8.0;

# Reconstruction via RED(ADMM)+SSA
admm_ssa,admm_it_ssa = red_admm!(fwd, adj, din,
                                 zero(dadj), λ, β, proj!, (dt,fmin,fmax,rank),
                                 ideal=d0,
                                 verbose=true,
                                 max_iter_o=K,
                                 max_iter_i1=10,
                                 max_iter_i2=1,
                                 tol=1e-6);



N=length(admm_it_ssa[:snr])
admm_qssa_snrx = zeros(N);
admm_qssa_snry = zeros(N);
admm_qssa_snrz = zeros(N);

for i in 1:N
    admm_qssa_snrx[i]=admm_it_ssa[:snr][i][1]
    admm_qssa_snry[i]=admm_it_ssa[:snr][i][2]
    admm_qssa_snrz[i]=admm_it_ssa[:snr][i][3]
end

## Plots
figname ="qrecon_quality_xyz_qssa"
figure(figname,figsize=(10,3))

subplot(131)
plot(1:L, pgd_qssa_snrx,label="PGD (QSSA) X")
plot(1:M, fp_qssa_snrx,label="REF-FP (QSSA) X")
plot(1:N, admm_qssa_snrx,label="RED-ADMM (QSSA) X")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(a)")
legend()


subplot(132)
plot(1:L, pgd_qssa_snry,label="PGD (QSSA) Y")
plot(1:M, fp_qssa_snry,label="REF-FP (QSSA) Y")
plot(1:N, admm_qssa_snry,label="RED-ADMM (QSSA) Y")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(b)")
legend()

subplot(133)
plot(1:L, pgd_qssa_snrz,label="PGD (QSSA) Z")
plot(1:M, fp_qssa_snrz,label="REF-FP (QSSA) Z")
plot(1:N, admm_qssa_snrz,label="RED-ADMM (QSSA) Z")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(c)")
legend()

tight_layout()

##########################################################################
function fk_thresh(IN::AbstractArray,sched::AbstractFloat)

    out = copy(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad
    tmp = PadOp(out; nin=n, npad=npad, flag="fwd")
    
    # fft: I do not have a fft! version for the qft but I do have a
    # fft overwrite with side="left" and μ = 1/√3(i,j,k) default.
    tmp .= fft(tmp)

    # threshold
    threshold!(tmp,sched)
    
    # Truncate
    out .= PadOp( ifft(tmp); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

# Define a patching operator
function proj!(state, (sched))
    # output allocation
    out = copy(state.x)
    
    # get iteration:
    it = state.it;
    
    # fk_thresh
    out .= fk_thresh(out,sched[it])
    
    # Return
    return out
end

fwd(x) = interp_ks3d(x,htt,h,3,10,"fwd")
adj(x) = interp_ks3d(x,htt,h,3,10,"adj")

dadj = adj(din);

# Step-size selection
α = 0.1;

# threshold scheduler
Pi, Pf, K = 99.9, 1, 101
sched = _sched(dadj,K,Pi,Pf,"exp") ./ 10 ;

# Reconstruction via PGD+FKT ≡ POCS
out_fkt,it_fkt = pgdls!(fwd, adj, din,
                        zero(dadj), proj!, (sched);
                        ideal=d0,
                        α = α, verbose=true,
                        maxIter=K, ε=1e-6);

pgd_qfkt_snrx = zeros(K);
pgd_qfkt_snry = zeros(K);
pgd_qfkt_snrz = zeros(K);

for i in 1:K
    pgd_qfkt_snrx[i]=it_fkt[:snr][i][1]
    pgd_qfkt_snry[i]=it_fkt[:snr][i][2]
    pgd_qfkt_snrz[i]=it_fkt[:snr][i][3]
end

# Reg param
λ = 5.0;

# Reconstruction via RED(FP)+FKT
red_fkt,red_it_fkt = red_fp!(fwd, adj, din,
                             zero(dadj), λ, proj!, (sched),
                             ideal=d0,
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-6);

fp_qfkt_snrx = zeros(K);
fp_qfkt_snry = zeros(K);
fp_qfkt_snrz = zeros(K);

for i in 1:K
    fp_qfkt_snrx[i]=red_it_fkt[:snr][i][1]
    fp_qfkt_snry[i]=red_it_fkt[:snr][i][2]
    fp_qfkt_snrz[i]=red_it_fkt[:snr][i][3]
end


# Reg param
λ = 5.0;

# ADMM-param
β = 5.0;

# Reconstruction via RED(ADMM)+FKT
admm_fkt,admm_it_fkt = red_admm!(fwd, adj, din,
                                 zero(dadj), λ, β, proj!, (sched),
                                 ideal=d0,
                                 verbose=true,
                                 max_iter_o=K,
                                 max_iter_i1=10,
                                 max_iter_i2=1,
                                 ε=1e-6);

admm_qfkt_snrx = zeros(K);
admm_qfkt_snry = zeros(K);
admm_qfkt_snrz = zeros(K);

for i in 1:K
    admm_qfkt_snrx[i]=admm_it_fkt[:snr][i][1]
    admm_qfkt_snry[i]=admm_it_fkt[:snr][i][2]
    admm_qfkt_snrz[i]=admm_it_fkt[:snr][i][3]
end


## Plots
figname ="qrecon_quality_xyz"
figure(figname,figsize=(10,3))

subplot(131)
plot(1:K, pgd_qfkt_snrx,label="PGD (QFT) X")
plot(1:K, fp_qfkt_snrx,label="REF-FP (QFT) X")
plot(1:K, admm_qfkt_snrx,label="RED-ADMM (QFT) X")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(a)")
legend()


subplot(132)
plot(1:K, pgd_qfkt_snry,label="PGD (QFT) Y")
plot(1:K, fp_qfkt_snry,label="REF-FP (QFT) Y")
plot(1:K, admm_qfkt_snry,label="RED-ADMM (QFT) Y")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(b)")
legend()


subplot(133)
plot(1:K, pgd_qfkt_snrz,label="PGD (QFT) Z")
plot(1:K, fp_qfkt_snrz,label="REF-FP (QFT) Z")
plot(1:K, admm_qfkt_snrz,label="RED-ADMM (QFT) Z")

xlabel("Iterations"*L" (j) ")
ylabel("Quality (dB)")
title("(c)")
legend()

tight_layout()
