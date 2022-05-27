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
