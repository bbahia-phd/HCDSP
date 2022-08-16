pwd()

using Distributed
addprocs(5)

@everywhere dev_dir=joinpath(homedir(),"projects")
@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()

@everywhere using Revise
@everywhere using HCDSP
@everywhere using SeisMain, SeisPlot, PyPlot
@everywhere using LinearAlgebra, FFTW, DelimitedFiles
@everywhere using Statistics

# Set path
work =  joinpath(homedir(),"Desktop/data/seismic_data/valhall/");
seis_files=joinpath(work,"seis");

# data sizes for 100x100 binning
elt=Float32; nt=2000; n1=110; n2=70; dt=elt(0.004f0);

# read seis file after binning
d,h,e = SeisRead(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_100x100.seis"));

# Clip to see later reflections a little... NB: This is no good
val = quantile(abs.(vec(d)), 1)
#d = elt.(clamp.(d, -val, val));

# sampling opetaror
S = SamplingOp(d);

# These are already new coords so I am rewritting stuff
sx = SeisMain.ExtractHeader(h,"sx"); sy = SeisMain.ExtractHeader(h,"sy");
gx = SeisMain.ExtractHeader(h,"gx"); gy = SeisMain.ExtractHeader(h,"gy");

# minimum to set up origin of reg grid
sx_max = maximum(sx);  sx_min = minimum(sx);
sy_max = maximum(sy);  sy_min = minimum(sy);
dsx = dsy = 100;

# Estimate travel-times for shift
TT = zeros(elt,n1,n2);
δt = 0.4; # for taper in seconds

vw = 1480; hw = 70; hw2 = hw*hw;

for itr in eachindex(sx)
    
    tsx = sx[itr];  tsy = sy[itr];

    tgx = gx[itr]; tgy = gy[itr];

    isx = floor(Int,(tsx - sx_min)/dsx)+1;
    isy = floor(Int,(tsy - sy_min)/dsy)+1;

    dsq = sqrt((tsx-tgx)^2+hw2)
    TT[isx,isy] = dsq/vw;
end

dT = zero(TT);
dT .= TT #.- minimum(TT);

# conventional record length
rec_length = dt*(nt-1);

# conv acqusition time in seconds (68hrs)
t_conv = rec_length*n1*n2

# number of saillines for each boat
nboat = n2 ÷ 2;

# boat1 fires 14 sailines regularly
tb1 = collect(0.0 : rec_length : rec_length*(n1*nboat)-1)

# boat2 fires 15 saillines randomly with respect to boat1
tmp = round.(Int64, (tb1 .+ 4.0 .* rand(n1*nboat))./dt);
tb2 = tmp.*dt;

# check i-th sample
i = 10
[tb1[i] tb2[i]]

# A grid to map sources and firing times
grid = Array{Tuple{Int},2}(undef,n1,n2);
grid = [(ix,iy)  for ix in 1:n1, iy in 1:n2];

# Split the grid into four boats
boat1 = Array{Tuple{Int},2}(undef,n1,nboat); boat1 = copy(grid[:,1:nboat]);
boat2 = Array{Tuple{Int},2}(undef,n1,nboat); boat2 = copy(grid[:,nboat+1:2nboat]);

# the firing times associated with the positions
nsx = zeros(Int,n1*n2); nsy = zeros(Int,n1*n2); tau = zeros(Float32,n1*n2);

for i in 1:n1*nboat
    tau[(1+(i-1)*2):(2+(i-1)*2)]=[tb1[i]; tb2[i]]
    nsx[(1+(i-1)*2):(2+(i-1)*2)]=[boat1[i][1]; boat2[i][1]]
    nsy[(1+(i-1)*2):(2+(i-1)*2)]=[boat1[i][2]; boat2[i][2]]
end

# blended time is 17 hours (add rec_length to listen last shot)
t_blend = maximum(tau) + rec_length;

# blending factor (approx 2)
beta = t_conv / t_blend;

# define blending forward and adjoint operators
PARAM = (nt = nt,
         nx = n1,
         ny = n2,
         dt = dt,
         tau = tau,
         sx = nsx,
         sy = nsy);

# TODO: Add the sampling operator here
fwd(x) = SeisBlendOp(S .* x,PARAM,"fwd");
adj(x) = S .* SeisBlendOp(x,PARAM,"adj");

# inverse crime: blend
b = fwd(d);

# pseudo-deblend
db = adj(b);

# Patching
psize = nextpow.(2,(500,10,10));
polap = (10,10,10);
 smin = (1,1,1);
 smax = (nt,n1,n2);

# for schedule
dpatch,pid = fwdPatchOp(db,psize,polap,smin,smax);

# threshold parameters
@everywhere Pi, Pf, N = 99.9, 0.01, 101;

# threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"abma") ./ 10;

# initial guess for all methods
d0 = zero(d);

####################################
###### First stage inversion #######

# PGD step-size (< 1/β ≈ 0.5)
α = elt(0.125);

# tolerance
ε = elt(1e-16);

# Deblending by inversion with non-robust denoiser
tmp,tmp_it = pgdls!(fwd, adj, b, d0,
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = d, α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results: pgd_fkt represents the inverted u1 in Li et al (conoco)
pgd_fkt = tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

#################################
##### start processing step #####

# Shift for best SSA
OUTP = copy(pgd_fkt);
OUT = ndshift(OUTP,dT,dt);

patches,pid = fwdPatchOp(OUT,psize,polap,smin,smax);

# define SSA function
@everywhere begin
    dt = 0.004f0;
    fmin = 0.0f0;
    fmax = 60.0f0;
    rank = 3;
    fssa(δ) = fx_process(δ,dt,fmin,fmax,fast_ssa_lanc,(rank))
end

# SSA patches
patches .= pmap(fssa,patches);

# Rewrite the solution
OUT = adjPatchOp(patches, pid, psize, polap, smin, smax);

# Band-passing
OUT .= elt.(reshape(SeisBPFilter(reshape(OUT,nt,:),dt,0,10,50,60),(nt,n1,n2)));

# Undo shift after SSA
OUT2 = ndshift(OUT,-dT, dt);

#######################
###### Debiasing ######
n = size(OUT);

# (Full-data) Fourier coefficients
α̂ = fft(S .* OUT);

# get mask 
M = zeros(eltype(α̂),size(α̂));
mα = median(abs.(α̂)) ./ 2;
gtmean(x) = x > mα
M[findall(gtmean,abs.(α̂))] .= one(eltype(M));

# Overall fwd and adj operators with transform
FwdOp(s) = S .* real(ifft(M .* s));
AdjOp(s) = M .* fft(S .* s)  ./ prod(n);

fwdDb(x) = fwd(FwdOp(x));
adjDb(x) = AdjOp(adj(x));

# initial model
x = M .* α̂;

# model blended data from new coefficients
bb = fwdDb(x); #FwdOp(x);

# residual
r = b .- bb; misfit = real(dot(r,r));

# gradient
g = adjDb(r); #AdjOp(r);
gprod = real(dot(g,g));

# conj grad
p = copy(g);

# max iter for debias step
max_iter = 100;

# flag
verbose = true;

for i in 1:max_iter
    q = fwdDb(p); #FwdOp(p);

    γ = gprod / (real(dot(q,q))+1e-10);

    # model and residual update
    x .+= γ .* p
    r .-= γ .* q
    misfit = real(dot(r,r));

    # grad
    g = adjDb(r); #AdjOp(r);

    gprod_new = real(dot(g,g));
    β = gprod_new / (gprod + 1e-10);
    gprod = copy(gprod_new)

    # conj grad
    p .= g .+ β .* p

    verbose ? println("Iteration $i misfit = $(misfit) and gradient = $(gprod)") : nothing
end

OUTP = FwdOp(x);

####################################
###### Second stage inversion ######

# under some processing, pgd_fkt has an estimate of the direct arrival
# so p1 is the blended direct arrival
p1 = fwd(S .* OUTP); 

# and u1 is the residual from the observed blended data
u1 = b .- p1;

# pseudo-deblend
du1 = adj(u1);

# for schedule
dpatch,pid = fwdPatchOp(du1,psize,polap,smin,smax);

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.01, 201, 10;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"exp") ./ 10;

### Second stage inversion is all based on u1 now ###

# Deblending by inversion with robust denoiser
tmp,tmp_it = pgdls!(fwd, adj, u1, d0,
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = (d .- OUT2), α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results
pgd_rfkt = tmp;
it_pgd_rfkt_snr = tmp_it[:snr];
it_pgd_rfkt_mis = tmp_it[:misfit];
