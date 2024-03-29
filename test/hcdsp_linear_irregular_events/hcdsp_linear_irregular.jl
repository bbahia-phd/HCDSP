
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

# kill me, kill me now, but this will give us access to wigb plotting 
using MATLAB
mat"addpath('./../../../QSEIS/src/Plotting/')"

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
mat"wigb($(di[:,1,:]),1,$x,$t)"

# decimate
perc1=50;

# indexes
indx0=1:nx*ny                           # all indexes
indx1=HCDSP.decimate_indexes(dzz,perc1) # missing indexes
indx2=setdiff(indx0,indx1)              # observed indexes

# (random) decimating regular data (data regular observed = dro)
dro = reshape(dr,nt,:);
dro[:,indx1] .= zero(eltype(dro));

# (random) decimating irregular data (data irregular observerd = dio)
dio = reshape(di,nt,:);
dio[:,indx1] .= zero(eltype(dio));

diy = copy(reshape(diy,nt,:));
diy[:,indx1] .= zero(eltype(diy));

dix = copy(reshape(dix,nt,:));
dix[:,indx1] .= zero(eltype(dix));
 
# Sampling operator (no need)
T = SamplingOp(dio);

# extract the real irregular without binning
dio2 = dio[:,indx2];
 
# how is this binning dio2: dio2 is a bunch of (observed) traces together in matrix form (albeit 3D)
# You are placing dio2 into regular bins in the dio array. This would be the binning step, since dio was gridded accordingly.
dio[:,indx2] .= dio2;

# The 3D versions of the regular and irregular decimated data
dr3d = reshape(dro,nt,nx,ny);
di3d = reshape(dio,nt,nx,ny);

di3dx = reshape(dix,nt,nx,ny);
di3dy = reshape(diy,nt,nx,ny);
di3dz = reshape(diz,nt,nx,ny);

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
    out .= fx_process(out,dt,fmin,fmax,fast_ssa_lanc,(rank))
    
    # Return
    return out
end

fwd(x) = interp_ks3d(x,htt,h,3,10,"fwd")
adj(x) = interp_ks3d(x,htt,h,3,10,"adj")

# Step-size selection
α = 0.1;

# Number iterations
K = 101;

# f-x process
dt=0.004; fmin=0; fmax=80; rank=10;

#d0 = copy(quaternion(dzx,dzy,dzz));

# Reconstruction via PGD+SSA (I-MSSA)
dadj = adj(dix2);
x_out_ssa,x_it_ssa = pgdls!(fwd, adj, dix2,
                          zero(dadj), proj!, (dt,fmin,fmax,rank),
                          ideal=dzx,
                          α = α, verbose=true,
                          maxIter=K, ε=1e-6);   

dadj = adj(diy2);
y_out_ssa,y_it_ssa = pgdls!(fwd, adj, diy2,
                          zero(dadj), proj!, (dt,fmin,fmax,rank),
                          ideal=dzy,
                          α = α, verbose=true,
                          maxIter=K, ε=1e-6);   

dadj = adj(diz2);
z_out_ssa,z_it_ssa = pgdls!(fwd, adj, diz2,
                          zero(dadj), proj!, (dt,fmin,fmax,rank),
                          ideal=dzz,
                          α = α, verbose=true,
                          maxIter=K, ε=1e-6);   

##########################################################################
# Reg param
λ = 2.5;

# Reconstruction via RED(FP)+SSA
red_ssa,red_it_ssa = red_fp!(fwd, adj, din,
                             zero(dadj), λ, proj!, (dt,fmin,fmax,rank),
                             ideal=d0,
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-6);

red_ssa,red_it_ssa = red_fp!(fwd, adj, din,
                             zero(dadj), λ, proj!, (dt,fmin,fmax,rank),
                             ideal=d0,
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-6);

red_ssa,red_it_ssa = red_fp!(fwd, adj, din,
                             zero(dadj), λ, proj!, (dt,fmin,fmax,rank),
                             ideal=d0,
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-6);

##########################################################################
# Reg param
λ = 2.5;

# ADMM-param
β = 2.5;

# Reconstruction via RED(ADMM)+SSA
admm_ssa,admm_it_ssa = red_admm!(fwd, adj, din,
                                 zero(dadj), λ, β, proj!, (dt,fmin,fmax,rank),
                                 ideal=d0,
                                 verbose=true,
                                 max_iter_o=K,
                                 max_iter_i1=10,
                                 max_iter_i2=1,
                                 tol=1e-6);

##########################################################################
function fk_thresh(IN::AbstractArray,sched::AbstractFloat)

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



#=
function fx_pgd_recon(in::AbstractArray{T},fwd::Function,adj::Function,Pi::Real,Pf::Real,K::Int) where {T} 

    # takes as input the frequency slice and prepare the fwd and adj operators
    dadj = adj(in); # this is a frequency slice nx × ny

    nx,ny = size(dadj);

    # threshold scheduler
    sched = _sched(dadj,K,Pi,Pf,"exp");

    # Step-size selection
    α = 0.4;

#=
    out = copy(dadj);

    for i in 1:K
        out .=  fk_thresh(out,sched[i]);

        yy1 = fwd(out);
        yy2 = adj(yy1);

        out .= dadj .+ out .- yy2

    end    
=#
    # Deblending by inversion    
    out,_  = pgdls!(fwd, adj, in,
                    zero(dadj), proj!, (sched);
                    ideal=zero(dadj),
                    α = α, verbose=false,
                    maxIter=K, tol=1e-6);   

    return out
end

fmin=0; fmax=80; 

out = HCDSP.fx_irregular_process(din,dt,fmin,fmax,htt,h,fx_pgd_recon,(fwd,adj,Pi,Pf,K)...);
=#

#=
# The dot product test of the interpolator
function interp_dot_prod_test()
    fwd(x) = interp_ks3d(x,htt,h,3,10,"fwd")
    adj(x) = interp_ks3d(x,htt,h,3,10,"adj")

    in1 = randn(size(dr3d)); # this is regular (undersampled) data
    in2 = randn(size(dio2)); # this is irregular (undersampled) data

    out1 = fwd(in1) # this turns in1 into irregular data
    out2 = adj(in2) # this turns in2 into regular data

    dot1 = dot(vec(in1),vec(out2))
    dot2 = dot(vec(in2),vec(out1))

    isapprox(dot1,dot2) ? println("Pass") : println("Fail")
end
   
interp_dot_prod_test()
=#

#=
# jitter: formal definition
# The basic idea of jittered undersampling is to regularly decimate the interpolation grid, thus generating a coarse grid, and subsequently perturb the coarse-grid sample points on the fine grid.

# As for random undersampling according to a discrete uniform distribution, where each location is equally likely to be sampled, a discrete uniform distribution for the peturbation around the coarse-grid points is considered.

# γ repreesnts the undersampling factor, here taken to be odd, γ = 1, 3, 5, ...
γx,γy = 5,5;

# SizeN of the interpolation gird is assumed to be a multiple of γ so that the number of acquired points n = N/γ is an integer.
nx = 10;
Nx = γx*nx;

ny = 10;
Ny = γy*ny;

# Jittered-sampled data points given by
# y[i] = f_0[j] for i = 1,...,n
# and j = 0.5*(1-γ) + γ * i + ϵ_i
# where ϵ_i is an iid (assume distribution) integer between [-floor(ξ-1)/2,floor(ξ-1)/2].
# The jitter parameter, 0 <= ξ <= γ, relates to the size of the perturbation around the corase-grid.
ξx,ξy = 5,5;

jx = zeros(Int,nx)
jy = zeros(Int,ny)

for i in 1:nx
    lx = -floor(0.5*(ξx-1))
    ux =  floor(0.5*(ξx-1))
    ϵxi = rand(lx:ux)

    ly = -floor(0.5*(ξy-1))
    uy =  floor(0.5*(ξy-1))
    ϵyi = rand(ly:uy)

    jx[i] = 0.5*(1-γx) + γx * i + ϵxi
    jy[i] = 0.5*(1-γy) + γy * i + ϵyi
end

# The case above is just "jittered sampling". If you set ξ=0, there is no jitter, and the undersampling scheme is just regular undersampling.
# The work of Hennenfent and Herrmann show the impact of ξ.
# Optimally-jittered undersampling sets the jitter parameter equal to the undersampling parameter: ξ = γ. Such a configuration is optimal because it creates the most favourable contitions for recovery with a localized transform.
=#
