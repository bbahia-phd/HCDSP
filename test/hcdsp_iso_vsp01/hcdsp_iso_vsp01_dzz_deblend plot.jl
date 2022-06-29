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
dx_path = joinpath(dir_path,"iso_vsp01_zx.seis");
dy_path = joinpath(dir_path,"iso_vsp01_zy.seis");
dz_path = joinpath(dir_path,"iso_vsp01_zz.seis");

fname = joinpath(dir_path,"all_iter_hist_fkt.h5")

fid = h5open(fname,"r")

it_pgd_fkt_snr  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_snr   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_snr = Vector{Vector{Float32}}(undef,nr);

it_pgd_fkt_mis  = Vector{Vector{Float32}}(undef,nr);
it_fp_fkt_mis   = Vector{Vector{Float32}}(undef,nr);
it_admm_fkt_mis = Vector{Vector{Float32}}(undef,nr);
for i in 1:nr
    it_admm_fkt_snr[i] = fid["admm"]["snr"]["$i"] |> read
    it_fp_fkt_snr[i]   = fid["fp"]["snr"]["$i"]   |> read
    it_pgd_fkt_snr[i]  = fid["pgd"]["snr"]["$i"]  |> read

    it_admm_fkt_mis[i] = fid["admm"]["mis"]["$i"] |> read
    it_fp_fkt_mis[i]   = fid["fp"]["mis"]["$i"]   |> read
    it_pgd_fkt_mis[i]  = fid["pgd"]["mis"]["$i"]  |> read
end

close(fid)


