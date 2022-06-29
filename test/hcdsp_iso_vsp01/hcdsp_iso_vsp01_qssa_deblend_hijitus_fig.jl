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

# read pseudo-deblended data
dx_path = joinpath(dir_path,"pseudo/iso_vsp01_zx_pseudo.seis");
dy_path = joinpath(dir_path,"pseudo/iso_vsp01_zy_pseudo.seis");
dz_path = joinpath(dir_path,"pseudo/iso_vsp01_zz_pseudo.seis");
dbx,_,_ = SeisRead(dx_path);
dby,_,_ = SeisRead(dy_path);
dbz,_,_ = SeisRead(dz_path);

# read ideal input data
dx_path = joinpath(dir_path,"inputs/iso_vsp01_zx.seis");
dy_path = joinpath(dir_path,"inputs/iso_vsp01_zy.seis");
dz_path = joinpath(dir_path,"inputs/iso_vsp01_zz.seis");
dzx,hzx,ext = SeisRead(dx_path);
dzy,hzy,ext = SeisRead(dy_path);
dzz,hzz,ext = SeisRead(dz_path);

# get geometry
sx = SeisMain.ExtractHeader(hzz,"sx");
sy = SeisMain.ExtractHeader(hzz,"sy");
gl = SeisMain.ExtractHeader(hzz,"gelev");

dt = 0.012;
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

ozx = zeros(eltype(dzx),nt,ns,nsline,nr);
ozy = zeros(eltype(dzy),nt,ns,nsline,nr);
ozz = zeros(eltype(dzz),nt,ns,nsline,nr);

for i in eachindex(hzz)

    # i-th trace header
    h = hzz[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    ozx[:,sxi,syi,gel] += dzx[:,i]
    ozy[:,sxi,syi,gel] += dzy[:,i]
    ozz[:,sxi,syi,gel] += dzz[:,i]
end

dzx = []; dzy = []; dzz = [];

include("./HCDSP/test/hcdsp_iso_vsp01/get_figs_func.jl")