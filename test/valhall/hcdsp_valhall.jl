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
ntb = 6207044;         # number of time samples in blended data
dt  = elt(0.004);      # sampling interval
nt  = 2000;            # time samles
nx  = 107;             # sources in x
ny  = 115;             # sources in y
 n  = (nt,nx,ny)       # tuple for inputs

data_path = joinpath(dev_dir,"HCDSP/test/valhall")

# read
tmp = read_write(joinpath(data_path,"shooting_times_638-738.bin"),"r",n=(3*nx*ny,1), T=elt);
tmp = reshape(tmp,(nx*ny,3));
tau = tmp[:,1]; sx = round.(Int,tmp[:,2]); sy = round.(Int,tmp[:,3]);

# Sampling operator
T = read_write(joinpath(data_path,"sampling_638-738.bin"),"r",n=(nt,54,58),T=elt);   
T = reshape(T,(nt,54,58));

# ground-truth data
d = read_write(joinpath(data_path,"sailline_638-738.bin"),"r",n=(nt,54*58));   
d = reshape(d,(nt,nx,ny));

# blended data
b = read_write(joinpath(data_path,"blended_638-738.bin"),"r",n=(ntb,1),T=elt);       

# Pseudo-deblend
PARAM = (nt = nt,     # time samples
         nx = nx,     # sources in x
         ny = ny,     # sources in y
         dt = dt,     # sampling in time
         tau = tau,   # firing times
         sx = sx,     # ordered list of shots x
         sy = sy);    # ordered list of shots y

db = SeisBlendOp(b,PARAM,"adj");

# Patching
psize = nextpow.(2,(500,20,20));
polap = (50,50,50);
 smin = (1,1,1);
 smax = (nt,nx,ny);

d2,pid = fwdPatchOp(db,psize,polap,smin,smax);

# Get threshold schedule
Pi, Pf, K = 99.0,0.01,101

sched = thresh_sched(d2,K,Pi,Pf,"exp") ./ 10;

figure("Schedule",figsize=(3,2.5))
plot(sched); gcf()

# define a projection operator
function proj!(state, (psize, polap, smin, smax, sched))
    
    # extract current iteration
    it = state.it;
   
    # apply patching on input
    patches,pid = fwdPatchOp(state.x, psize, polap, smin, smax);
    
    # apply fk_thresh to all arrays in tmp using the threshold at iteration it
    patches = pmap(fk_thresh, patches, repeat([sched[it]], length(patches)));
    
    # rewrite the solution
    state.x = adjPatchOp(patches, pid, psize, polap, smin, smax);
   
    # return
    #return state 
end

@everywhere function fk_thresh(IN::AbstractArray, t::Real)

    out = similar(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)
  
    # pad
    tmp = complex.(PadOp(IN; nin=n, npad=npad, flag="fwd"))
    
    # Analysis
    fft!(tmp)
    
   # Threshold
    threshold!(tmp,t)
   
    # Synthesis 
    ifft!(tmp);
    
    # truncate
    out = PadOp(real(tmp); nin=n, npad=npad, flag="adj")
    
    return out
end

# Pseudo-deblend
bFwd(x) = SeisBlendOp(x, PARAM, "fwd");
bAdj(x) = SeisBlendOp(x, PARAM, "adj");

α = elt(0.5);
ε = elt(1e-6);

# zeros of same type
d0 = zero(db);

@time x_pgd,it = pgdls!(bFwd, bAdj, b, d0,
                        proj!, (psize,polap,smin,smax,sched);
                        ideal = d0,  maxIter=K, α = α, ε = ε);



