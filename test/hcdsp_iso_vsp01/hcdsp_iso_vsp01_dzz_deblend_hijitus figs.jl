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

include(joinpath(dev_dir,"dev/deblending/src/deblend_module.jl"))
import Main.deblend_module: BlendOp

using PyPlot
using SeisMain, SeisPlot

using HDF5

# data dir home
dir_path  = "/media/bbahia/DATA/seismic_data/iso_vsp01"

# read input data
dzz,hzz,ext_zz = SeisRead(joinpath(dir_path,"iso_vsp01_ozz.seis"));
db,_,_         = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pseudo.seis"));


pgd_fkt ,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pgd_fkt.seis"));
fp_fkt  ,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_fp_fkt.seis"));
admm_fkt,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_admm_fkt.seis"));

# get geometry
sx = SeisMain.ExtractHeader(hzz,"sx");
sy = SeisMain.ExtractHeader(hzz,"sy");
gl = SeisMain.ExtractHeader(hzz,"gelev");

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;

nsline = 205; # number of source lines
ns = 205;     # number of sources within lines
nt = 218;     # number of time samples
nr = 31;      # number of receivers down-well
drz = 16.6 ;  # source and receiver spacing (m )

dt = ext_zz.d1;
taxis = collect(0:nt-1).*dt;

# inputs
get_fig_1()
get_fig_4()

# residuals fkt
r1 = dzz .- pgd_fkt;
r2 = dzz .- fp_fkt;
r3 = dzz .- admm_fkt;

include("./HCDSP/test/hcdsp_iso_vsp01/get_figs_func.jl")

get_fig_2()
get_fig_3(h=12,w=10)

# ssa
pgd_fkt ,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pgd_ssa.seis"));
fp_fkt  ,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_fp_ssa.seis"));
admm_fkt,_,_   = SeisRead(joinpath(dir_path,"iso_vsp01_zz_admm_ssa.seis"));

# residuals ssa
r1 = dzz .- pgd_fkt;
r2 = dzz .- fp_fkt;
r3 = dzz .- admm_fkt;

get_fig_2()
get_fig_3(h=12,w=10)
