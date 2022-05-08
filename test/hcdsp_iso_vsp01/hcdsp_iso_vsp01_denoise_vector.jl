pwd()

using Distributed
addprocs(19)

@everywhere dev_dir="/dev/Breno_GOM/projects/";

cd(dev_dir)
#cd(joinpath(homedir(),"projects"))

pwd()

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP/"))
@everywhere Pkg.status()

@everywhere using Revise

@everywhere using FFTW
@everywhere using PyPlot
@everywhere using LinearAlgebra
@everywhere using SeisMain, SeisPlot
@everywhere using HCDSP#, IterativeMethods

# data dir home
data_path = "/dev/Breno_GOM/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zx_crg1350.seis");
dzy,hzy,ext_zy = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zy_crg1350.seis");
dzz,hzz,ext_zz = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zz_crg1350.seis");

dt      = Float64(ext_zz.d1);
nt      = Int64(ext_zz.n1);
ntr     = ext_zz.n2;
rz_init = 1350.0; # m
rz_end  = 1850.01; # m
drz     = 16.667 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr      = length(rz_axis)
ns      = 205; # number of sources within lines
nsline  = 205; # number of source lines

@assert ns*nsline == ntr == size(dzz,2)

dzx = reshape(dzx,(nt,ns,nsline));
dzy = reshape(dzy,(nt,ns,nsline));
dzz = reshape(dzz,(nt,ns,nsline));

# Add noise
snrx,snry,snrz=0.8,1.0,1.2;
dnx = SeisAddNoise(dzx, snrx, db=true);
dny = SeisAddNoise(dzy, snry, db=true);
dnz = SeisAddNoise(dzz, snrz, db=true);

# Missing traces
perc = 50;

# Temporary Quaternion
Qc = quaternion(dzx,dzy,dzz);
Qt = quaternion(dnx,dny,dnz);

# Decimate quaternion: all have same traces missing
Qt .= decimate_traces(Qt,perc);

@everywhere begin
    dt = 0.012012;
    fmin = 0;
    fmax = 50;
    psize = (64,32,32);
    polap = (50,50,50);
    smin = (1,1,1);
    smax = (217,205,205);
end

dx = imagi.(Qt);
dy = imagj.(Qt);
dz = imagk.(Qt);

# apply patching on input
pc,_        = fwdPatchOp(Qc, psize, polap, smin, smax);
patches,pid = fwdPatchOp(Qt, psize, polap, smin, smax);

# Define operator to act on a frequency slice d
@everywhere imp_ssa(d,k) = HCDSP.imputation_op(d, fast_qssa_lanc, (k); iter=50)

# "rank"-selection
@everywhere k = 15;

# define fx ssa function
@everywhere fssa(δ) = fx_process(δ,dt,fmin,fmax,imp_ssa,(k))

# time single-patch
@time pp = fssa(patches[10]);

# SSA all patches
pqout = similar(patches);
@time pqout .= pmap(fssa,patches);    

# patch difference
difp = pc .- pqout;

# unpatch the solution
xout = adjPatchOp(imagi.(pqout), pid, psize, polap, smin, smax);
yout = adjPatchOp(imagj.(pqout), pid, psize, polap, smin, smax);
zout = adjPatchOp(imagk.(pqout), pid, psize, polap, smin, smax);

# check scale difference
α = dot(dzz,zout) / dot(zout,zout)
dif = dzz - α .* zout;

#n=100;
#SeisPlotTX([dzz[:,:,n] dz[:,:,n] out[:,:,n] dif[:,:,n]],
#            wbox=8, hbox=4, cmap="gray");
#gcf()
