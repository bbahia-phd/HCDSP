pwd()

using Distributed
addprocs(5)

cd(joinpath(homedir(),"projects/HCDSP/test/hyper_events"))

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

@everywhere using Revise
@everywhere using FFTW
@everywhere using HCDSP, LinearAlgebra, Random

using SeisMain, SeisPlot, SeisProcessing
using PyPlot

# same as collect(ot:dt:tmax) or collect(range(start::T,stop::T,step::T=))
ot = 0.0; dt = 0.004; nt = 301; tmax = ot + (nt-1)*dt
taxis = range(ot,tmax,length=nt) |> collect; 

# set spatial dimension
ox = -1000.0; dx = 20; nx = 101; xmax = ox + (nx-1)*dx;
xaxis = range(ox,xmax,length=nx) |> collect;

# zero-offset travel-times
tau = [0.2, 0.6, 0.9]; 

# rms velocities
vel = [1500, 2000, 3000];

# apex-shifts (for dipping layers)
apex = zero(vel);

# events amplitudes
amp = [5, -0.1, 0.1];

# wavelet
wav_flag = "ricker";

# central freq
f0 = [20.0];

d = SeisHypEvents(f0 = f0,
                wavelet =  wav_flag,
                amp  = amp,
                apex = apex,
                vel  = vel,
                tau  = tau,
                ot  = ot,
                ox  = ox,
                nx  = nx,
                nt  = nt );
SeisPlotTX(d,pclip=95); gcf()

# seed
Random.seed!(1992)

# blending factor
β = 0.5;

# total time 
t_blended = tmax*β*nx;

tau = sort(rand(nx)*t_blended)

# set random source positions
sx = 1:nx |> collect |> shuffle!
sy = ones(Int,size(sx));

# Pseudo-deblend
PARAM = (nt = nt,     # time samples
         nx = nx,     # sources in x
         ny =  1,     # sources in y
         dt = dt,     # sampling in time
         tau = tau,   # firing times
         sx = sx,     # ordered list of shots x
         sy = sy);    # ordered list of shots y

# Pseudo-deblend
bFwd(x) = SeisBlendOp(x, PARAM, "fwd")[:,1];
bAdj(x) = SeisBlendOp(x, PARAM, "adj")[:,:,1];

# I will work with a single receiver
b = bFwd(d);
db = bAdj(b);
dc = copy(d);

pgd_fkt  = similar(db);
fp_fkt   = similar(db);
admm_fkt = similar(db);

# Threshold schedule
@everywhere Pi, Pf, K = 99.9, 0.1, 101

# Params for patching
psize = nextpow.(2,(20,20));
polap = (10,10);
smin = (1,1);
smax = (nt,nx);

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

# Select a CRG
tmpz = dc;

# Inverse crime: model observed data
b  = bFwd(tmpz);

# Pseudo-deblended
db = bAdj(b);

# Patching
dpatch,_ = fwdPatchOp(db,psize,polap,smin,smax);

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,K,Pi,Pf,"exp");

####################################

# PGD step-size (< 1/β ≈ 0.5)
α = Float64(0.4);

# Deblending by inversion
tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = dc, α=α,
                    verbose=true,
                    maxIter=K,
                    ε=Float64(1e-16));

# Store inversion results
pgd_fkt = tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

####################################

# RED reg param
λ = Float64(0.5);

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
λ = Float64(0.5); γ = Float64(0.5);

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

SeisPlotTX(tmp,pclip=99,vmin=-a,vmax=a,cmap="gray")




SeisPlotTX([db dc],pclip=99); gcf()