cd(joinpath(homedir(),"projects"))
pwd()

using Distributed
addprocs(3)

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
@everywhere Pkg.status()

@everywhere using Revise

@everywhere using FFTW
@everywhere using PyPlot
@everywhere using SeisMain, SeisPlot
@everywhere using HCDSP, IterativeMethods

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz_crg1350.seis"));

dt = Float64(ext_zz.d1);
nt = Int64(ext_zz.n1);
ntr = ext_zz.n2;
rz_init = 1350.0; # m
rz_end  = 1850.01; # m
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
snrx,snry,snrz=0.8,1.0,1.2;
dnx = SeisAddNoise(dzx, snrx, db=true);
dny = SeisAddNoise(dzy, snry, db=true);
dnz = SeisAddNoise(dzz, snrz, db=true);

# Missing traces
perc = 40;

# Temporary Quaternion
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
pc,_        = fwdPatchOp(dzz, psize, polap, smin, smax);
patches,pid = fwdPatchOp(dz, psize, polap, smin, smax);

# Define operator to act on a frequency slice d
@everywhere imp_ssa(d,k) = HCDSP.imputation_op(d, fast_ssa_lanc, (k); iter=5, α=0.4)

# "rank"-selection
@everywhere k = 8;

# define fx ssa function
@everywhere fssa(δ) = fx_process(δ,dt,fmin,fmax,imp_ssa,(k))

@time pp = fssa(patches[10]);

# fk_thresh all patches
pout = similar(patches);
@time pout .= pmap(fssa,patches);    

# rewrite the solution
out = adjPatchOp(pout, pid, psize, polap, smin, smax);
dif = dzz - out;

n=100;
SeisPlotTX(
    [dzz[:,:,n] dz[:,:,n] out[:,:,n] dif[:,:,n]],wbox=8, hbox=4, cmap="gray");
gcf()