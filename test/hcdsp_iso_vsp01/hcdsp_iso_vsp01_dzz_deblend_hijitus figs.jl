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
dir_path = "/media/bbahia/DATA/seismic_data/iso_vsp01"

# read input data
dzz,hzz,ext_zz = SeisRead(joinpath(dir_path,"inputs/iso_vsp01_ozz.seis"));
dbz,_,_        = SeisRead(joinpath(dir_path,"pseudo/iso_vsp01_zz_pseudo.seis"));

pgd_fkt ,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pgd_fkt.seis"));
fp_fkt  ,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_fp_fkt.seis"));
admm_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_admm_fkt.seis"));

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

include("./HCDSP/test/hcdsp_iso_vsp01/get_figs_func.jl")

# residuals fkt
r1 = 5 .* (dzz .- fp_fkt);
r2 = 5 .* (dzz .- pgd_fkt);
r3 = 5 .* (dzz .- admm_fkt);

get_fig_2()
get_fig_3(h=12,w=10)

# ssa
pgd_fkt ,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_pgd_ssa.seis"));
fp_fkt  ,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_fp_ssa.seis"));
admm_fkt,_,_ = SeisRead(joinpath(dir_path,"iso_vsp01_zz_admm_ssa.seis"));

# residuals ssa
r1 = dzz .- pgd_fkt;
r2 = dzz .- fp_fkt;
r3 = dzz .- admm_fkt;

get_fig_2()
get_fig_3(h=12,w=10)

# inputs
get_fig_1()
get_fig_4()


fname = joinpath(dir_path,"fkt/crg15/all_iter_hist_fkt.h5")

fid = h5open(fname,"r")
nr = 1;

it_pgd_fkt_snr  = Vector{Float32}(undef,K);
it_fp_fkt_snr   = Vector{Float32}(undef,K);
it_admm_fkt_snr = Vector{Float32}(undef,K);

it_pgd_fkt_mis  = Vector{Float32}(undef,K);
it_fp_fkt_mis   = Vector{Float32}(undef,K);
it_admm_fkt_mis = Vector{Float32}(undef,K);

#=
it_fp_fkt_snr   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Float32}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float32}}(undef,nr);
=#

for i in 1:nr
    it_admm_fkt_snr[i] = fid["admm"]["snr"]["$i"] |> read
    it_fp_fkt_snr[i]   = fid["fp"]["snr"]["$i"]   |> read
    it_pgd_fkt_snr[i]  = fid["pgd"]["snr"]["$i"]  |> read

    it_admm_fkt_mis[i] = fid["admm"]["mis"]["$i"] |> read
    it_fp_fkt_mis[i]   = fid["fp"]["mis"]["$i"]   |> read
    it_pgd_fkt_mis[i]  = fid["pgd"]["mis"]["$i"]  |> read
end

close(fid)


