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

# sampling and samples
elt = Float32;
ntb = 1620510;#6207044;   # number of time samples in blended data
dt  = elt(0.004);         # sampling interval
nt  = 2000;               # time samles
nx  = 54;                 # sources in x
ny  = 58;                 # sources in y
 n  = (nt,nx,ny)          # tuple for inputs

data_path = joinpath(homedir(),"Desktop/data/seismic_data/valhall/bin");

# read
tmp = read_write(joinpath(data_path,"shooting_times_638-738.bin"),"r",n=(3*nx*ny,1), T=elt);
tmp = reshape(tmp,(nx*ny,3));
tau = tmp[:,1]; sx = round.(Int,tmp[:,2]); sy = round.(Int,tmp[:,3]);

# Sampling operator
T = read_write(joinpath(data_path,"sampling_638-738.bin"),"r",n=(nt,54,58),T=elt);   
T = reshape(T,(nt,nx,ny));

# ground-truth data
d = read_write(joinpath(data_path,"sailline_638-738.bin"),"r",n=(nt,54,58),T=elt);   
d = reshape(d,(nt,nx,ny));

# blended data
b = read_write(joinpath(data_path,"blended_638-738.bin"),"r",n=(ntb,1),T=elt);       

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
    αi = Float32(0.5);

    # iht tolerance
    εi = Float32(1e-6);    

    # Robust Iterative Hard Thresholding
    tmp,_ = riht!(FwdOp, AdjOp, out, zeros(ComplexF64,npad),
                  sched, update_weights, p;
                  #α = αi,
                  maxIter=K,
                  verbose=false)

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
    #new_sched = range(sched[1], stop=sched[it], length=K)
    
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

# Pseudo-deblend
PARAM = (nt = nt,     # time samples
         nx = nx,     # sources in x
         ny = ny,     # sources in y
         dt = dt,     # sampling in time
         tau = tau,   # firing times
         sx = sx,     # ordered list of shots x
         sy = sy);    # ordered list of shots y

# Pseudo-deblended
db = SeisBlendOp(b,PARAM,"adj");

# Patching
psize = nextpow.(2,(200,10,10));
polap = (50,50,50);
 smin = (1,1,1);
 smax = (nt,nx,ny);

# for schedule
dpatch,pid = fwdPatchOp(db,psize,polap,smin,smax);

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.1, 201, 10;

# Threshold scheduler
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"abma") ./ 10;

figure("Schedule",figsize=(3,2.5))
plot(sched); gcf()

# Pseudo-deblend
bFwd(x) = SeisBlendOp(x, PARAM, "fwd");
bAdj(x) = SeisBlendOp(x, PARAM, "adj");

# p-vals for robust thresholding
p = [1.6, 1.7,1.8,1.9,2.0];
nintervals = length(p);
Ni = div(N,nintervals);

pvals = zeros(Float64,N);
c, cc = 0,0;
for j in 1:nintervals
    global c += 1;
    for i in 1:Ni
        global cc += 1;
        pvals[cc] = p[c]
    end
end
pvals[cc+1]=2.0;

# initial guess for all methods
d0 = zero(d);

####################################
# PGD step-size (< 1/β ≈ 0.5)
α = elt(0.25);

# tolerance
ε = elt(1e-16);

# Deblending by inversion with non-robust denoiser
tmp,tmp_it = pgdls!(bFwd, bAdj, b, d0,
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = d, α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results
pgd_fkt = tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

# Deblending by inversion with robust denoiser
tmp,tmp_it = pgdls!(bFwd, bAdj, b, d0,
                    rproj!, (psize,polap,smin,smax,sched,pvals);
                    ideal = d, α = α,
                    verbose=true,
                    maxIter=N,
                    ε=ε);

# Store inversion results
pgd_rfkt = tmp;
it_pgd_rfkt_snr = tmp_it[:snr];
it_pgd_rfkt_mis = tmp_it[:misfit];

####################################
# RED reg param
λ = Float64(0.1);

# Deblending by inversion with robust denoiser
tmp,tmp_it = red_fp!(bFwd, bAdj, b, zero(db), λ,
                     rproj!, (psize,polap,smin,smax,sched,pvals);
                     ideal = tmpz,
                     verbose=true,
                     max_iter_o=N,
                     max_iter_i=10,
                     ε=Float32(1e-16));

# Store inversion results
fp_rfkt = tmp;
it_fp_rfkt_snr = tmp_it[:snr];
it_fp_rfkt_mis = tmp_it[:misfit];
