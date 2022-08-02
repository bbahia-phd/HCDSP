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
amp = [10, -1, 0.1];

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
b = bFwd(d); db = bAdj(b); dc = copy(d);

# Threshold parameters
@everywhere Pi, Pf, N, K = 99.9, 0.1, 201, 100;

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

# Params for patching
psize = nextpow.(2,(20,20));
polap = (50,50);
smin = (1,1);
smax = (nt,nx);

########################################################################
@everywhere function update_weights(W,r,pvals)#; γ=2.0)
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
        #γ2=γ^2
        # @inbounds for i in eachindex(W)
        #     W[i] = γ2 / ( γ2 +  r[i]^2.0);
        # end
        @inbounds for i in eachindex(W)
            W[i] = 1 / (abs(r[i])^(2.0-p) +  ε)
        end        
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

    # Pad
    # tmp = complex.( PadOp(IN; nin=n, npad=npad, flag="fwd") )
    
    # Operator handles
    #FwdOp(s) = ifft(s) ;
    #AdjOp(s) =  fft(s) ./ prod(npad);

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

    # smooth traces
    #out .= smooth_traces(out);

    # mute prearrival
    #out[1:400,:] .= 0.0
    
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
    
    # smooth traces
    #out .= smooth_traces(out);

    # mute prearrival
    #out[1:400,:] .= 0.0
    
    # projection
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
sched = HCDSP.thresh_sched(dpatch,N,Pi,Pf,"abma") ./ 10;

####################################
# PGD step-size (< 1/β ≈ 0.5)
α = Float64(0.1);

# Deblending by inversion
tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                    proj!, (psize,polap,smin,smax,sched);
                    ideal = dc, α=α,
                    verbose=true,
                    maxIter=N,
                    ε=Float64(1e-16));

# Store inversion results
pgd_fkt = tmp;
it_pgd_fkt_snr = tmp_it[:snr];
it_pgd_fkt_mis = tmp_it[:misfit];

# Deblending by inversion 
tmp,tmp_it = pgdls!(bFwd, bAdj, b, zero(db),
                    rproj!, (psize,polap,smin,smax,sched,pvals);
                    ideal = dc, α=α,
                    verbose=true,
                    maxIter=N,
                    ε=Float64(1e-16));

                    # Store inversion results
pgd_rfkt = tmp;
it_pgd_rfkt_snr = tmp_it[:snr];
it_pgd_rfkt_mis = tmp_it[:misfit];

####################################
# RED reg param
λ = Float64(0.1);

# Deblending by inversion    
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

####################################
j=50;
tmp = [d db pgd_fkt pgd_rfkt fp_fkt];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=90,vmin=-a,vmax=a,cmap="gray")

j=50;
tmp = [d db 5 .* (d .- pgd_fkt) 5 .* (d .- pgd_rfkt)];
a = maximum(tmp[:])*0.2;

SeisPlotTX(tmp,pclip=99,vmin=-a,vmax=a,cmap="gray")

SeisPlotTX([db dc],pclip=99); gcf()

####################################

figname="rfkt_snr_plot";
close("all")
figure(figname,figsize=(8,8) .* 0.5)

plot(1:N,it_pgd_fkt_snr, label="PGD (FKT)",     color="black")
plot(1:N,it_pgd_rfkt_snr, label="PGD (RFKT)",    color="red")
plot(1:N,it_fp_rfkt_snr, label="RED-FP (RFKT)", color="green")
title("(a)")
xlabel("Iteration number (k)",fontsize=15)
ylabel(L"R_z = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_2 }{\parallel \hat{\bf d}_k - {\bf d}^{o} \parallel^2_2} \right)}",fontsize=15)
ylim(0,80)
xticks([40;80;120;160;200])
legend()

tight_layout()
gcf()

savefig(joinpath("/home/bbahia/projects/files",figname))
