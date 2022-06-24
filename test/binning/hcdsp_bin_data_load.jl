import Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using Revise
using PyPlot
using HCDSP

using MAT, MATLAB
mat"addpath('./../../../QSEIS/src/Plotting/')"

using SeisMain, SeisPlot, SeisProcessing

# Doing this so SeisMain doesn't send things away
ENV["DATAPATH"]=""

# Check mat file
file = matopen("data_patch.mat")
varnames = keys(file)

# Read header
head = read(file,"H");

dt  = head["dt"][1];

sx = head["sx"][1,:] |> Vector{Float64};
sy = head["sy"][1,:] |> Vector{Float64};

gx = head["gx"][1,:] |> Vector{Float64};
gy = head["gy"][1,:] |> Vector{Float64};

dout = Float32.(read(file,"D")[800:1400,:]);
dout = SeisBPFilter(dout,dt / 1e6 ,1,3,50,65);

nt,ntrace = size(dout);
t = collect((0:nt-1) .* dt) ./ 1e6;

x0 = min(minimum(sx),minimum(gx)) + 10
y0 = min(minimum(sy),minimum(gy)) + 10

# translate origin
sx .-= x0
sy .-= y0
gx .-= x0
gy .-= y0

# survey azimuth
sr_az = 39π/180;
c = cos(sr_az); s = sin(sr_az);

# rotated source position
sxr = sx .* c .- sy .* s;
syr = sx .* s .+ sy .* c;

# rotated receiver pos
gxr = gx .* c .- gy .* s;
gyr = gx .* s .+ gy .* c;

# for mid points
mx = similar(sxr)
my = similar(syr)

# for offsets
hx = similar(sxr)
hy = similar(syr)

# for azimuth
h = similar(sxr)
az = similar(syr)

# some fields are left blank (0.0)
for n in 1:ntrace
   
    # mid points x and y
    mx[n] = (sxr[n]+gxr[n])/2;
    my[n] = (syr[n]+gyr[n])/2;

    # offsets
    hx[n] = (sxr[n]-gxr[n]);
    hy[n] = (syr[n]-gyr[n]);
    h[n]  = sqrt(hx[n]^2+hy[n]^2);

    az[n] = (180/π) * atan(hy[n],hx[n])
    if az[n] < 0
        az[n] += 360.0
    end

end

# define a regular grid
mx_min = minimum(mx); mx_max = maximum(mx);
my_min = minimum(my); my_max = maximum(my);
h_max = maximum(h); h_min = minimum(h);
az_min = 0; az_max =  2*180;

dmx = 20;
nx = floor(Int,(mx_max-mx_min)/dmx) + 1;

dmy = 20;
ny = floor(Int,(my_max-my_min)/dmy) + 1;

dh = 400;
nh = floor(Int,(h_max - h_min)/dh) + 1;

daz = 45;
naz = floor(Int,(az_max-az_min)/daz)+1;

# sampling operator
T = zeros(nx,ny,nh,naz); 

# bin counter
count = zeros(nx,ny,nh,naz);

# data bin
dbin = zeros(eltype(dout),nt,nx,ny,nh,naz);

# Initialize regular grids
mx_grid = [];
my_grid = [];
h_grid  = [];
az_grid = [];

for n in 1:ntrace
    i = floor(Int,(mx[n]-mx_min)/dmx)+1;
    j = floor(Int,(my[n]-my_min)/dmy)+1;
    k = floor(Int,(h[n]-h_min)/dh)+1;
    l = floor(Int,(az[n]-az_min)/daz)+1;


    T[i,j,k,l] = 1;
    count[i,j,k,l] += 1;

    dbin[:,i,j,k,l] .+= dout[:,n]

    append!(mx_grid, mx_min + (i-1) * dmx)
    append!(my_grid, my_min + (j-1) * dmy)
    append!(h_grid,  h_min  + (k-1) * dh)
    append!(az_grid, az_min + (l-1) * daz)
end
 
cindex = CartesianIndices((1:nx,1:ny,1:nh,1:naz));
tot = 0;
tot_bin =  prod((nx,ny,nh,naz));
for i in eachindex(T)
    if T[i] == 1
        tot += 1;
        dbin[:,cindex[i]] ./= count[i]
    end
end

trc_perc = round(tot/tot_bin*100,digits=2)
println("$tot out of $tot_bin ($trc_perc %) alive traces")

# Define a patching operator
function proj!(state, (dt,fmin,fmax,rank))
    # output allocation
    out = copy(state.x)
    
    # get iteration:
    it = state.it;
    
    # fk thresh
    out .= fx_process(out,dt,fmin,fmax,fast_ssa_lanc,(rank))
    
    # Return
    return out
end

# Overwrite the sampling operator
T = SamplingOp(dbin);

# set fwd/adj pair (notice adjoint as identity)
fwd(x) = T .* x;
adj(x) = x;

# Step-size selection
α = 0.1;

# Number iterations
K = 201;

# f-x process
dt  = head["dt"][1] / 1e6; fmin=0; fmax=65; rank=4;

# Reconstruction via PGD+SSA (I-MSSA)
out_ssa,it_ssa = pgdls!(fwd, adj, dbin, zero(dbin),
                        proj!, (dt,fmin,fmax,rank),
                        ideal=zero(dbin),
                        α = α, verbose=true,
                        maxIter=K, ε=1e-16);

# Reg param
λ = 8;

# Reconstruction via RED(FP)+SSA
red_ssa,red_it_ssa = red_fp!(fwd, adj, dbin, zero(dbin), λ,
                             proj!, (dt,fmin,fmax,rank),
                             ideal=zero(dbin),
                             verbose=true,
                             max_iter_o=K,
                             max_iter_i=10,
                             ε=1e-16);
