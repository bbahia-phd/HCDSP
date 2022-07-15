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

using MAT,MATLAB
mat"addpath('/home/bbahia/projects/QSEIS/src/Plotting/')"

using HDF5


# data dir home
dir_path  = "/media/bbahia/DATA/seismic_data/iso_vsp01/crg1350"

dzx,hzx,ezx = SeisRead(joinpath(dir_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzx,ezx = SeisRead(joinpath(dir_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzx,ezx = SeisRead(joinpath(dir_path,"iso_vsp01_zz_crg1350.seis"));

nzx,hzx,ezx = SeisRead(joinpath(dir_path,"noisy_iso_vsp01_zx_crg1350.seis"));
nzy,hzx,ezx = SeisRead(joinpath(dir_path,"noisy_iso_vsp01_zy_crg1350.seis"));
nzz,hzx,ezx = SeisRead(joinpath(dir_path,"noisy_iso_vsp01_zz_crg1350.seis"));

rzx,hzx,ezx = SeisRead(joinpath(dir_path,"ssa_iso_vsp01_zx_crg1350.seis"));
rzy,hzx,ezx = SeisRead(joinpath(dir_path,"ssa_iso_vsp01_zy_crg1350.seis"));
rzz,hzx,ezx = SeisRead(joinpath(dir_path,"ssa_iso_vsp01_zz_crg1350.seis"));

qzx,hzx,ezx = SeisRead(joinpath(dir_path,"qssa_iso_vsp01_zx_crg1350.seis"));
qzy,hzx,ezx = SeisRead(joinpath(dir_path,"qssa_iso_vsp01_zy_crg1350.seis"));
qzz,hzx,ezx = SeisRead(joinpath(dir_path,"qssa_iso_vsp01_zz_crg1350.seis"));

azx,hzx,ezx = SeisRead(joinpath(dir_path,"aqssa_iso_vsp01_zx_crg1350.seis"));
azy,hzx,ezx = SeisRead(joinpath(dir_path,"aqssa_iso_vsp01_zy_crg1350.seis"));
azz,hzx,ezx = SeisRead(joinpath(dir_path,"aqssa_iso_vsp01_zz_crg1350.seis"));

n = nt,n1,n2 = 217,205,205;
drz = 16.6;

dzx = permutedims(reshape(dzx,n),(1,3,2));
dzy = permutedims(reshape(dzy,n),(1,3,2));
dzz = permutedims(reshape(dzz,n),(1,3,2));

nzx = permutedims(reshape(nzx,n),(1,3,2));
nzy = permutedims(reshape(nzy,n),(1,3,2));
nzz = permutedims(reshape(nzz,n),(1,3,2));

rzx = permutedims(reshape(rzx,n),(1,3,2));
rzy = permutedims(reshape(rzy,n),(1,3,2));
rzz = permutedims(reshape(rzz,n),(1,3,2));

qzx = permutedims(reshape(qzx,n),(1,3,2));
qzy = permutedims(reshape(qzy,n),(1,3,2));
qzz = permutedims(reshape(qzz,n),(1,3,2));

azx = permutedims(reshape(azx,n),(1,3,2));
azy = permutedims(reshape(azy,n),(1,3,2));
azz = permutedims(reshape(azz,n),(1,3,2));

# get geometry
sx = SeisMain.ExtractHeader(hzx,"sx");
sy = SeisMain.ExtractHeader(hzx,"sy");
gl = SeisMain.ExtractHeader(hzx,"gelev");

sx_n,sy_n = sx |> copy,sy |> copy;

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;

for i in eachindex(hzx)

    # i-th trace header
    h = hzx[i]

    sxi = floor(Int,(h.sx - sx_min)/drz)+1
    syi = floor(Int,(h.sy - sy_min)/drz)+1
#    gel = floor(Int,(-h.gelev + gl_min)/drz)+1

    trace = nzx[:,sxi,syi];
    if norm(trace,2)^2 < 1e-8
        sx_n[i] = NaN
        sx_n[i] = NaN
    end    
   
end
