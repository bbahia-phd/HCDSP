pwd()

dev_dir = "/home/bbahia/projects";
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
dir_path  = "/media/bbahia/DATA/seismic_data/iso_vsp01/"

# read ideal input data
dz_path = joinpath(dir_path,"inputs/iso_vsp01_zx.seis");
dzx,hzz,ext = SeisRead(dz_path);
# read ideal input data
dz_path = joinpath(dir_path,"inputs/iso_vsp01_zy.seis");
dzy,hzz,ext = SeisRead(dz_path);
# read ideal input data
#dz_path = joinpath(dir_path,"inputs/iso_vsp01_zz.seis");
dz_path = joinpath(dir_path,"inputs/iso_vsp01_ozz.seis");
dzz,hzz,ext = SeisRead(dz_path);

# get geometry
sx = SeisMain.ExtractHeader(hzz,"sx");
sy = SeisMain.ExtractHeader(hzz,"sy");
gl = SeisMain.ExtractHeader(hzz,"gelev");

dt = 0.012f0;

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

# read pgd-qssa results
x_pgd_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zx_pgd_qfkt.seis"));
y_pgd_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zy_pgd_qfkt.seis"));
z_pgd_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zz_pgd_qfkt.seis"));

# read fp-qssa results
x_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zx_fp_qfkt.seis"));
y_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zy_fp_qfkt.seis"));
z_fp_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zz_fp_qfkt.seis"));

# read admm-qssa results
x_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zx_admm_qfkt.seis"));
y_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zy_admm_qfkt.seis"));
z_admm_fkt,_,_ = SeisRead(joinpath(dir_path,"qfkt/iso_vsp01_zz_admm_qfkt.seis"));


#=

r1 = dzz .- z_pgd_fkt;
r2 = dzz .- z_fp_fkt;
r3 = dzz .- z_admm_fkt;

=#


for i in eachindex(hzz)

    # i-th trace header
    h = hzz[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1
    
    x_pgd_fkt[:,sxi,syi,gel] .= 1 .* (dzx[:,i] .- x_pgd_fkt[:,sxi,syi,gel])
    y_pgd_fkt[:,sxi,syi,gel] .= 1 .* (dzy[:,i] .- y_pgd_fkt[:,sxi,syi,gel])
    z_pgd_fkt[:,sxi,syi,gel] .= 1 .* (dzz[:,i] .- z_pgd_fkt[:,sxi,syi,gel])

    x_fp_fkt[:,sxi,syi,gel] .= 1 .* (dzx[:,i] .- x_fp_fkt[:,sxi,syi,gel])
    y_fp_fkt[:,sxi,syi,gel] .= 1 .* (dzy[:,i] .- y_fp_fkt[:,sxi,syi,gel])
    z_fp_fkt[:,sxi,syi,gel] .= 1 .* (dzz[:,i] .- z_fp_fkt[:,sxi,syi,gel])

    x_admm_fkt[:,sxi,syi,gel] .= 1 .* (dzx[:,i] .- x_admm_fkt[:,sxi,syi,gel])
    y_admm_fkt[:,sxi,syi,gel] .= 1 .* (dzy[:,i] .- y_admm_fkt[:,sxi,syi,gel])
    z_admm_fkt[:,sxi,syi,gel] .= 1 .* (dzz[:,i] .- z_admm_fkt[:,sxi,syi,gel])
   
end

dzx = []; dzy = []; dzz=[];
