cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using PyPlot
using SeisMain, SeisPlot
using HCDSP, IterativeMethods

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz_crg1350.seis"));

dt = ext_zz.d1;
nt = ext_zz.n1;

ntr = ext_zz.n2;

rz_init = 1350.0f0; # m
rz_end  = 1850.01f0; # m
drz = 16.667 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr = length(rz_axis)

nsline = 205; # number of source lines
ns = 205; # number of sources within lines

@assert ns*nsline == ntr == size(dzz,2)

dzx = reshape(dzx,(nt,ns,nsline));
dzy = reshape(dzy,(nt,ns,nsline));
dzz = reshape(dzz,(nt,ns,nsline));

# Add noise
dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
dny = SeisAddNoise(dzy, -2.0, db=true, L=3);
dnz = SeisAddNoise(dzz, -2.0, db=true, L=3);

# Temporary Quaternion
Qt = quaternion(dnx,dny,dnz);

# Missing traces
Qt .= decimate_traces(Qt,perc);

fmin = 0; fmax = 80;
psize = (256,32,32); polap = (50,50,50);
smin = (1,1,1); smax = (nt,ns,nsline);

# apply patching on input
patches,pid = fwdPatchOp(Qt, psize, polap, smin, smax);

# Define operator to act on a frequency slice d
imp_ssa(d,k) = HCDSP.imputation_op(d,SVDSSAOp,(k);iter=10)

# define fx ssa function
fssa(δ) = fx_process(δ,dt,fmin,fmax,imp_ssa,(rank[it]))

# fk_thresh all patches
patches .= pmap(fssa,patches);    

# rewrite the solution
out = adjPatchOp(patches, pid, psize, polap, smin, smax);

# # Component-wise processing
# Xo = fx_process(imagi.(Qt),dt,fmin,fmax,imp_ssa,(k))
# Yo = fx_process(imagj.(Qt),dt,fmin,fmax,imp_ssa,(k))
# Zo = fx_process(imagk.(Qt),dt,fmin,fmax,imp_ssa,(k))

# # Call fx_process with Q imputation
# Qo = fx_process(Qt,dt,fmin,fmax,imp_ssa,(2k))

