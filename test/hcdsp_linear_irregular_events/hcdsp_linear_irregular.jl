
import Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using Revise
using LazyGrids
using HCDSP
using SeisPlot, SeisMain

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

# Ok. This is how Rongzhi is doing:
nx,ny,nt=30,30,200;
dx,dy,dt=20,20,0.004;

# midpoints
mx = collect(0:nx-1).*dx;
my = collect(0:ny-1).*dy;

(gmx,gmy) = ndgrid(mx,my)

jx = rand(size(gmx)...) .* 2 .- 1
jy = rand(size(gmy)...) .* 2 .- 1

jmx = gmx .+ jx; jmx[1,:] .= abs.(jmx[1,:]);
jmy = gmy .+ jy; jmy[:,1] .= abs.(jmy[:,1]);

px = [1/4000, 1/9000, -1/4000];
py = [1/4000, -1/4000, -1/6000];

τ  = [0.1,0.3,0.6];
amp = [1,-1,1];

f0 = 20.0;

d = HCDSP.SeisLinearIrregularEvents();
