pwd()

using Distributed
addprocs(30)

@everywhere dev_dir=joinpath(homedir(),"projects")
@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()

@everywhere using Revise
@everywhere using HCDSP
@everywhere using SeisMain, SeisPlot, SeisProcessing, PyPlot
@everywhere using LinearAlgebra, FFTW, DelimitedFiles
@everywhere using Statistics,Random
@everywhere include("./ndshift.jl")

# fixed-seed rng
rng = MersenneTwister(1992);

# Set path
work =  joinpath(homedir(),"projects/valhall");

# data sizes for 100x100 binning
elt=Float32; nt=2000; n1=110; n2=70; dt=elt(0.004f0);

# read time shifts here
dT = read_write(joinpath(homedir(),"projects/HCDSP/test/valhall/bin/travel_times.bin"),"r",n=(n1,n2),T=elt);
dT = reshape(dT,(n1,n2));

# read seis file after binning
d,h,e = SeisRead(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_100x100.seis"));
d = d[1:1000,1:50,1:50];
nt,n1,n2 = size(d);

# Clip to see later reflections a little... NB: This is no good
val = quantile(abs.(vec(d)), 0.9f0)
d  .= elt.(clamp.(d, -val, val));

# should I normalize?
cind = CartesianIndices((n1,n2));
for k in cind
    amax = maximum([abs(d[i,k]) for i in 1:nt])
    if amax == 0
	amax= 0.0001
    end
    d[:,k] .= d[:,k] ./ amax
end

# number of saillines for each boat
nboat = n2 ÷ 2;

# conventional record length
rec_length = dt*(nt-1);

# conv acqusition time in seconds (68hrs)
t_conv = rec_length*n1

# blending factor
β = 0.9;

# total time 
t_blended = rec_length*β*n1;

tau = elt.(sort(rand(n1)*t_blended))

# set random source positions
sx = 1:n1*n2 |> collect |> shuffle!
sy = ones(Int,size(sx));

# Pseudo-deblend
PARAM = (nt = nt,     # time samples
         nx = n1,     # sources in x
         ny =  1,     # sources in y
         dt = dt,     # sampling in time
         tau = tau,   # firing times
         sx = sx,     # ordered list of shots x
         sy = sy);    # ordered list of shots y

# boat1 fires 14 sailines regularly
tb1 = collect(0.0 : rec_length : rec_length*(n1*nboat)-1);

# boat2 fires 15 saillines randomly with respect to boat1
tmp = round.(Int64, (tb1 .+ 4.0 .* rand(rng,n1*nboat))./dt);
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

# Sampling opetaror
S = SamplingOp(d);

# time gain
#tg,itg = get_gain(d,dt;a=1.0);
#tg  .= sqrt.(tg);
#itg .= sqrt.(itg);

# TODO: Add the sampling operator here
fwd(x) = SeisBlendOp(x, PARAM, "fwd");
adj(x) = SeisBlendOp(x, PARAM, "adj");

# inverse crime: blend
b = fwd(d);

# pseudo-deblend
db = adj(b);

# Patching
psize = nextpow.(2,(20,20,20));
polap = (50,50,50);
 smin = (1,1,1);
 smax = (nt,n1,n2);

# for schedule
dpatch,pid = fwdPatchOp(db,psize,polap,smin,smax);

# Threshold parameters
@everywhere Pi, Pf, N = 99, 0.01, 101;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"exp") ./ 20;

# initial guess for all methods
d0 = zero(d);

# tolerance
ε = elt(1e-16);

# RED trade-off
λ = elt(0.25);

# Deblending by inversion
tmp,tmp_it = red_fp!(fwd, adj, b, d0, λ,
                     proj!, (psize,polap,smin,smax,sched);
                     ideal = d,
                     verbose=true,
                     max_iter_o=N,
                     max_iter_i=50,
                     ε=ε);

# Store inversion results
rfp_fkt = tmp;
it_rfp_fkt_snr = tmp_it[:snr];
it_rfp_fkt_mis = tmp_it[:misfit];

# p-vals for robust thresholding
p = [1.8,1.9,2.0];
nintervals = length(p);
Ni = div(N,nintervals);

pvals = zeros(elt,N);
c,cc = 0,0;
for j in 1:nintervals
    global c += 1;
    for i in 1:Ni
        global cc += 1;
        pvals[cc] = p[c]
    end
end
pvals[cc+1] = 2.0f0;

# RED trade-off
λ = elt(0.25);

# Spectral refinement
@everywhere K = 5;

# Deblending by inversion with robust deinoisers
tmp,tmp_it = red_fp!(fwd, adj, b, d0, λ,
                     rproj!, (psize,polap,smin,smax,sched,pvals);
                     ideal = d,
                     verbose=true,
                     max_iter_o=N,
                     max_iter_i=5,
                     ε=ε);

# Store inversion results
rfp_rfkt = tmp;
it_rfp_rfkt_snr = tmp_it[:snr];
it_rfp_rfkt_mis = tmp_it[:misfit];

j = 20; pc=99;
to_plot = [d[:,:,j] db[:,:,j] rfp_fkt[:,:,j] rfp_rfkt[:,:,j]];

SeisPlotTX(to_plot,cmap="gray",
           vmin=-val,vmax=val,
           pclip=pc,name="tmp",
           wbox=12,
           hbox=5,
#           xticks=[55,165,275,385],
#           xticklabels=[L"{\bf d}",L"{\bf d}_b","FKT","RFKT"],
           oy=0.0,
           dy=dt,
           ylabel="Time (s)");
