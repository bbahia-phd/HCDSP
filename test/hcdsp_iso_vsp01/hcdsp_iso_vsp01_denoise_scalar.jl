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
rz_init = 1350.0;  # m
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

@everywhere begin
    
    # f-x params
    dt    = 0.012012;
    fmin  = 0;
    fmax  = 50;
    
    # patching params
    psize = (64,32,32);
    polap = (50,50,50);
    smin  = (1,1,1);
    smax  = (217,205,205);

    # rank-reduction params
    k     = 15;
    iter  = 50;

    # Define operator to act on a frequency slice d
    imp_ssa(d,k) = HCDSP.imputation_op(d,
                            fast_ssa_lanc,
                            (k); iter=iter)

    # define fx ssa function to act on time domian δ
    fssa(δ) = fx_process(δ,dt,fmin,fmax,imp_ssa,(k))
end

function pmap_fx_ssa(δ,δc)
    # apply patching on input
    cl_patch, pid = fwdPatchOp(δc, psize, polap, smin, smax);
    in_patch, _   = fwdPatchOp(δ, psize, polap, smin, smax);

    # SSA all patches
    out_patch  = similar(in_patch);
    out_patch .= pmap(fssa,in_patch);    

    # un-patch
    out_full = adjPatchOp(out_patch, pid, psize, polap, smin, smax);

    # scale out to clean
    α = dot(δc,out_full) / dot(out_full,out_full)
    
    # get differences
    diff_patch = cl_patch .- out_patch;
    diff_full  = δc .- α .* out_full;    

    return out_full,diff_full,out_patch,diff_patch,α,cl_patch,in_patch
end

# Add noise
snrx,snry,snrz=0.8,1.0,1.2;
dnx = SeisAddNoise(dzx, snrx, db=true,L=3);
dny = SeisAddNoise(dzy, snry, db=true,L=3);
dnz = SeisAddNoise(dzz, snrz, db=true,L=3);

# Decimate quaternion: all have same traces missing
perc = 50;
Qt = quaternion(dnx,dny,dnz);
Qt = decimate_traces(Qt,perc);

# split inputs
dx = imagi.(Qt);
dy = imagj.(Qt);
dz = imagk.(Qt);

x_full,x_diff_full,x_patch,x_diff_patch,x_α,x_cl_patch,x_in_patch = pmap_fx_ssa(dx, dzx)
y_full,y_diff_full,y_patch,y_diff_patch,y_α,y_cl_patch,y_in_patch = pmap_fx_ssa(dy, dzy)
z_full,z_diff_full,z_patch,z_diff_patch,z_α,z_cl_patch,z_in_patch = pmap_fx_ssa(dz, dzz)

close("all");clf()
j = 10; n=10;
SeisPlotTX([x_cl_patch[j][:,:,n] x_patch[j][:,:,n] x_diff_patch[j][:,:,n]],
            xcur=2,
            style="overlay",
            hbox=4, wbox=12);
gcf()

n=100;
SeisPlotTX([dzz[:,:,n] dz[:,:,n] out[:,:,n] dif[:,:,n]], wbox=8,  hbox=4, cmap="gray");
gcf()

n=120;
SeisPlotTX([dzz[n,:,:] dz[n,:,:] out[n,:,:] dif[n,:,:]], wbox=10, hbox=3, cmap="gray");
gcf()