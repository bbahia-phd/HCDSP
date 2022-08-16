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

# Set path
work =  joinpath(homedir(),"Desktop/data/seismic_data/valhall/");
bin_files = joinpath(work,"bin");
su_files=joinpath(work,"su");
seis_files=joinpath(work,"seis");
figs=joinpath(work,"figs");

file_name=readdir(su_files)
file_name=file_name[ findall( x -> occursin(".su",x), file_name ) ]

data_path=joinpath.(su_files,file_name[ findall( x -> occursin(".su",x), file_name ) ])
save_path=similar(data_path);

# data sizes for 100x100 binning
elt=Float32; nt=2000; n1=110; n2=70; dt=elt(0.004f0);

# read seis file after binning
d,h,e = SeisRead(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_100x100.seis"));

# sampling opetaror
S = SamplingOp(d);

# These are already new coords so I am rewritting stuff
sx = SeisMain.ExtractHeader(h,"sx"); sy = SeisMain.ExtractHeader(h,"sy");
gx = SeisMain.ExtractHeader(h,"gx"); gy = SeisMain.ExtractHeader(h,"gy");

# minimum to set up origin of reg grid
sx_max = maximum(sx);  sx_min = minimum(sx);
sy_max = maximum(sy);  sy_min = minimum(sy);

# Estimate travel-times for shift
TT = zeros(elt,n1,n2); TTup = zero(TT); TTdw = zero(TT);
δt = 0.4; # for taper in seconds

vw = 1500; hw = 70; hw2 = hw*hw;

for itr in eachindex(TT)

    tsx = sx[itr]; tsy = sy[itr];
    tgx = gx[itr]; tgy = gy[itr];

    dsq = sqrt((tsx-tgx)^2+(tsy-tgy)^2+hw2)
    TT[itr] = dsq/vw;
    TTup[itr] = TT[itr]-δt
    TTdw[itr] = TT[itr]+δt
end

dT = zero(TT);
dT .= TT .- minimum(TT);

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

########################################################################
@everywhere function update_weights(W,r,pvals,ii; γ=2.0)
    # The function riht! multiplies W element by element.
    ε = 1e-4;
   
    p = pvals

    if p == 1.0
        @inbounds for i in eachindex(W)
            W[i] = 1 / ( abs(r[i]) +  ε);
        end

    elseif p == 2.0
        @inbounds for i in eachindex(W)
            W[i] = one(r[i]);
        end
    else
        @inbounds for i in eachindex(W)
            W[i] = 1 / (abs(r[i])^(2.0-p) +  ε)
        end
        #γ2=γ^2
        # @inbounds for i in eachindex(W)
        #     W[i] = γ2 / ( γ2 +  r[i]^2.0);
        # end        
    end
end

########################################################################
# robust thresholding
@everywhere function fk_thresh(IN::AbstractArray,sched::AbstractArray, p)

    out = copy(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad & Crop
    fwdPad(x) = PadOp(x; nin=n, npad=npad, flag="fwd");
    adjPad(x) = PadOp(x; nin=n, npad=npad, flag="adj");
    
    # Overall fwd and adj operators with transform
    FwdOp(s) = adjPad(real(ifft(s)));
    AdjOp(s) = fft(fwdPad(s))  ./ prod(npad);    

    # iht step-size
    αi = Float32(0.1);

    # iht tolerance
    εi = Float32(1e-16);    

    # Robust Iterative Hard Thresholding
    tmp,_ = riht!(FwdOp, AdjOp, out, zeros(ComplexF32,npad),
                  sched, update_weights, p;
#                  α = αi,
                  maxIter=K,
                  verbose=false)

    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end


########################################################################
# define a patching operator
function rproj!(state, (psize, polap, smin, smax, sched, pvals))

    # output allocation
    out = similar(state.x);
    
    # get iteration:
    it = state.it;

    # p value
    p = pvals[it];

    # set internal sched for RIHT
    new_sched = _schedule(sched[1], sched[it], K, "exp")
    
    # apply patching on input
    patches,pid = fwdPatchOp(state.x, psize, polap, smin, smax);
    
    # fk_thresh all patches
    if p != 2.0
        patches .= pmap(fk_thresh,
                   patches,
                   repeat([new_sched], length(patches)),
                   repeat([p],length(patches)));
    else
        # define the fkt function with given thresholding                                  
        fkt(δ) = fk_thresh(δ,sched[it])
        patches .= pmap(fkt,patches);
    end
        
    # rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);
        
    # projection
    return out
end

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

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.01, 151, 10;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"exp") ./ 10;

#figure("Schedule",figsize=(3,2.5))
#plot(sched); gcf()

# initial guess for all methods
d0 = zero(d);

####################################
###### First stage inversion #######

# PGD step-size (< 1/β ≈ 0.5)
α = elt(0.125);

# tolerance
ε = elt(1e-16);

# p-vals for robust thresholding
p = [1.6,1.7,1.8,1.9,2.0];
nintervals = length(p);
Ni = div(N,nintervals);

pvals = zeros(elt,N);
c, cc = 0,0;
for j in 1:nintervals
    global c += 1;
    for i in 1:Ni
        global cc += 1;
        pvals[cc] = p[c]
    end
end
pvals[cc+1]=2.f0;

# Deblending by inversion with robust denoiser
tmp,tmp_it = pgdls!(fwd, adj, b, d0,
                    rproj!, (psize,polap,smin,smax,sched,pvals);
                    ideal = d, α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results: pgd_fkt represents the inverted u1 in Li et al (conoco)
pgd_rfkt = tmp;
it_pgd_rfkt_snr = tmp_it[:snr];
it_pgd_rfkt_mis = tmp_it[:misfit];
