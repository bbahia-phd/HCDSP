pwd()

dev_dir = "/dev/Breno_GOM/projects";
cd(dev_dir)
pwd()

using Pkg
Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()

using Revise
using LinearAlgebra
using FFTW
using HCDSP

using PyPlot
using SeisMain, SeisPlot

using HDF5

# data dir home
dir_path  = "/dev/Breno_GOM/projects/files/iso_vsp01/"

# read pgd-qssa results
x_pgd_fkt,hzz,ext = SeisRead(joinpath(dir_path,"iso_vsp01_zx_pgd_qssa.seis"));
y_pgd_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zy_pgd_qssa.seis"));
z_pgd_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pgd_qssa.seis"));

# read fp-qssa results
x_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zx_fp_qssa.seis"));
y_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zy_fp_qssa.seis"));
z_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_fp_qssa.seis"));

# read admm-qssa results
x_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zx_admm_qssa.seis"));
y_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zy_admm_qssa.seis"));
z_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_admm_qssa.seis"));

# residuals fkt
r1 = ozz .- z_pgd_fkt;
r2 = ozz .- z_fp_fkt;
r3 = ozz .- z_admm_fkt;

# get geometry
sx = SeisMain.ExtractHeader(hzz,"sx");
sy = SeisMain.ExtractHeader(hzz,"sy");
gl = SeisMain.ExtractHeader(hzz,"gelev");

nt = 218;     # number of time samples (int64)
ns = 205;     # number of sources within lines
nsline = 205; # number of source lines
nr = 31;      # number of receivers down-well
drz = 16.6 ;  # source and receiver spacing (m )

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;
