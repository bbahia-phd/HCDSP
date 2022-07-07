pwd()

dev_dir = joinpath(homedir(),"projects")

using Pkg
Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()


# Usage packages
using PyPlot
using LinearAlgebra
using FFTW
using DelimitedFiles

using SeisPlot,SeisMain
using HCDSP

# save dir
data_path = "/media/bbahia/DATA/seismic_data/missi/missi_blend_fold_2/"

# sampling and samples
ntb = 692037; nt = 1751; nr = 183; ns = 790; dt = Float32(4e-3)

# read ground truth
d = read_write(joinpath(data_path,"conventional_missi.bin"),"r";n=(nt,ns,nr),T=Float32);
d = reshape(d,(nt,ns,nr));

# read blended data
bfull = read_write(joinpath(data_path,"blended_missi.bin"),"r";n=(ntb,nr),T=Float32);
bfull = reshape(bfull,(ntb,nr));

# read blending parameters
tmp = read_write(joinpath(data_path,"shooting_times_missi.bin"),"r";n=(3*ns,1),T=Float32);
tmp = reshape(tmp,ns,3);
tau = tmp[:,1]; sx = Int.(tmp[:,2]); sy = Int.(tmp[:,3]);

# Pseudo-deblend
PARAM = (nt = nt,     # time samples
         nx = ns,     # sources in x
         ny =  1,     # sources in y
         dt = dt,     # sampling in time
         tau = tau,   # firing times
         sx = sx,     # ordered list of shots x
         sy = sy);    # ordered list of shots y

# Pseudo-deblend
bFwd(x) = SeisBlendOp(x, PARAM, "fwd")[:,1];
bAdj(x) = SeisBlendOp(x, PARAM, "adj")[:,:,1];

# I will work with a single receiver
b = bfull[:,50];

db2 = bAdj(b);

d_clean = copy(d[:,:,50]);

w=6;
h=12

figure("Results",figsize=(w,h))

subplot(2,1,1)
SeisPlotTX(d_clean,
           fignum="Results",
           cmap="gray",
           pclip=90,
           dy=dt,
           ylabel="Time (s)",
           xlabel="Shot index",
           title="(a) Clean",
           ticksize=10,
           labelsize=14,
           titlesize=14)

subplot(2,1,2)
SeisPlotTX(db2,
           fignum="Results",
           cmap="gray",
           pclip=90,
           dy=dt,
           ylabel="Time (s)",
           xlabel="Shot index",
           title="(b) Combed",    
           ticksize=10,
           labelsize=14,
           titlesize=14)

tight_layout();

savefig("./input_crg_bf2.pdf")
