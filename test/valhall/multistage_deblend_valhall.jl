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
elt=Float32; nt=2000; n1=220; n2=140; dt=elt(0.004f0);

# Read seis file after binning
d,h,e = SeisRead(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_50x50.seis"));

# its sampling opetaror
T = read_write(joinpath(work,"bin/binned_sampling_sailline_638-738_50x50.bin"), "r"; n = (n1,n2), T = elt);
T = reshape(T,(n1,n2));

# conventional record length
rec_length = dt*(nt-1);

# conv acqusition time in seconds (68hrs)
t_conv = rec_length*n1*n2

# number of saillines for each boat
nboat = n2 ÷ 4;

# boat1 fires 14 sailines regularly
tb1 = collect(0.0 : rec_length : rec_length*(n1*nboat-1))

# boat2 fires 15 saillines randomly with respect to boat1
tmp = round.(Int64, (tb1 .+ 4.0 .* rand(n1*nboat))./dt);
tb2 = tmp.*dt;

# boat3 fires 15 saillines randomly with respect to boat2
tmp = round.(Int64, (tb2 .+ 4.0 .* rand(n1*nboat))./dt);
tb3 = tmp.*dt;

# boat4 fires 15 saillines randomly with respect to boat3
tmp = round.(Int64, (tb3 .+ 4.0 .* rand(n1*nboat))./dt);
tb4 = tmp.*dt;

# check i-th sample
i = 10
[tb1[i] tb2[i] tb3[i] tb4[i]]

# A grid to map sources and firing times
grid = Array{Tuple{Int},2}(undef,n1,n2)
grid = [(ix,iy)  for ix in 1:n1, iy in 1:n2]

# Split the grid into four boats (boat2 has less shots)
boat1 = Array{Tuple{Int},2}(undef,n1,nboat); boat1 = copy(grid[:,1:nboat]);
boat2 = Array{Tuple{Int},2}(undef,n1,nboat); boat2 = reverse(reverse(copy(grid[:,3nboat+1:4nboat]),dims=1),dims=2);
boat3 = Array{Tuple{Int},2}(undef,n1,nboat); boat3 = copy(grid[:,nboat+1:2nboat]);
boat4 = Array{Tuple{Int},2}(undef,n1,nboat); boat4 = reverse(reverse(copy(grid[:,2nboat+1:3nboat]),dims=1),dims=2);

# the firing times associated with the positions
nsx = zeros(Int,n1*n2); nsy = zeros(Int,n1*n2); tau = zeros(Float32,n1*n2);

for i in 1:n1*nboat
    tau[(1+(i-1)*4):(4+(i-1)*4)]=[tb1[i]; tb2[i]; tb3[i]; tb4[i]]
    nsx[(1+(i-1)*4):(4+(i-1)*4)]=[boat1[i][1]; boat2[i][1]; boat3[i][1]; boat4[i][1]]
    nsy[(1+(i-1)*4):(4+(i-1)*4)]=[boat1[i][2]; boat2[i][2]; boat3[i][2]; boat4[i][2]]
end

# blended time is 17 hours (add rec_length to listen last shot)
t_blend = maximum(tau) + rec_length;

# blending factor (approx 4)
beta = t_conv / t_blend

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
    tmp,_ = riht!(FwdOp, AdjOp, out, zeros(ComplexF64,npad),
                  sched, update_weights, p;
                  α = αi,
                  maxIter=K,
                  verbose=true)

    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

########################################################################
# non robust thresholding
@everywhere function fk_thresh(IN::AbstractArray,sched::Real)

    out = similar(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad
    tmp = complex.(PadOp(IN; nin=n, npad=npad, flag="fwd"))
    
    # fft
    fft!(tmp)

    # threshold
    threshold!(tmp,sched)
    
    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

########################################################################
# define a patching operator
function proj!(state, (psize, polap, smin, smax, sched))

    # output allocation
    out = similar(state.x);
    
    # get iteration:
    it = state.it;

    # apply patching on input
    patches,pid = fwdPatchOp(state.x, psize, polap, smin, smax);
    
    # define the fkt function with given thresholding                                                                                  
    fkt(δ) = fk_thresh(δ,sched[it])

    patches .= pmap(fkt,patches);
        
    # rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);

    # projection
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
# should be fwd(x) = SeisBlendOp(T .* x, PARAM, "fwd")
fwd(x) = SeisBlendOp(x,PARAM,"fwd");
adj(x) = SeisBlendOp(x,PARAM,"adj");

# inverse crime: blend
b = fwd(d);

# pseudo-deblend
db = adj(b);

# Patching
psize = nextpow.(2,(200,10,10));
polap = (50,50,50);
 smin = (1,1,1);
 smax = (nt,n1,n2);

# for schedule
dpatch,pid = fwdPatchOp(db,psize,polap,smin,smax);

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.1, 201, 10;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"exp") ./ 5;

figure("Schedule",figsize=(3,2.5))
plot(sched); gcf()

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

# Store inversion results:
# pgd_fkt represents the inverted u1 in Li et al (conoco)
pgd_fkt = tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

#################################
##### start processing step #####

# Shift traces by TT
# Padding
IN = copy(pgd_fkt);
nin  = size(IN);
npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

# Allocation
OUT = zero(IN);
OUTF = zero(INF);

# Fourier in time
fft!(INF,1);

# Freq range
ω_range = freq_indexes(0, 1000, dt, npad[1]);

# Spatial indexes
indx = CartesianIndices( npad[2:end] );
dw = 2π/npad[1]/dt;

# Loop over freqs
@inbounds for iω in ω_range
    # shift
    w = (iω-1)*dw;
    shift = exp.(1im*w .* dT[indx])
    OUTF[iω,indx] .= INF[iω,indx] .* shift
end

# Symmetries and ifft
conj_symmetry!(OUTF)

# Truncate
OUT .= PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );

patches,pid = fwdPatchOp(OUT,psize,polap,smin,smax);

# define ssa function
@everywhere begin
    dt = 0.004f0;
    fmin = 0.f0;
    fmax = 150.f0;
    rank = 3;
end

@everywhere fssa(δ) = fx_process(δ,dt,fmin,fmax,fast_ssa_qr,(rank))

# fk_thresh all patches
patches .= pmap(fssa,patches);

# Rewrite the solution
OUTP = adjPatchOp(patches, pid, psize, polap, smin, smax);

##### Undo Shift traces by TT
IN = copy(OUTP);
nin  = size(IN);
npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

# Allocation
OUT2 = zero(IN);
OUTF2 = zero(INF);

# Fourier in time
fft!(INF,1);

# Freq range
ω_range = freq_indexes(0, 150, dt, npad[1]);

# Spatial indexes
indx = CartesianIndices( npad[2:end] );
dw = 2π/npad[1]/dt;

# Loop over freqs
@inbounds for iω in ω_range
    # shift
    w = (iω-1)*dw;
    shift = exp.(-1im*w .* dT[indx])
    OUTF2[iω,indx] .= INF[iω,indx] .* shift
end

# Symmetries and ifft
conj_symmetry!(OUTF2)

# Truncate
OUT2 .= PadOp( real( ifft!(OUTF2,1) ), nin = nin, npad = npad, flag="adj" );

####################################
###### Second stage inversion ######

# under some processing, pgd_fkt has an estimate of the direct arrival
# so p1 is the blended direct arrival
p1 = fwd(OUT2); 

# and u1 is the residual from the observed blended data
u1 = b .- p1;

# pseudo-deblend
du1 = adj(u1);

# for schedule
dpatch,pid = fwdPatchOp(du1,psize,polap,smin,smax);

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.0001, 201, 10;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"exp") ./ 10;

figure("Schedule_2",figsize=(3,2.5))
plot(sched); gcf()

### Second stage inversion is all based on u1 now ###

# p-vals for robust thresholding
p = [1.6, 1.7,1.8,1.9,2.0];
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
tmp,tmp_it = pgdls!(fwd, adj, u1, d0,
                    rproj!, (psize,polap,smin,smax,sched,pvals);
                    ideal = (d .- OUT2), α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results
pgd_rfkt = tmp;
it_pgd_rfkt_snr = tmp_it[:snr];
it_pgd_rfkt_mis = tmp_it[:misfit];
