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
dir_path  = "/media/bbahia/DATA/seismic_data/linear5d";

n  = nt,n1,n2,n3,n4 = 100,20,20,20,20;
dt = 0.004;

# ideal data

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zx.bin")
dzx = read_write(file,"r",n=n);
dzx = reshape(dzx,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zy.bin")
dzy = read_write(file,"r",n=n);
dzy = reshape(dzy,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zz.bin")
dzz = read_write(file,"r",n=n);
dzz = reshape(dzz,n);

# noisy data
file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zx.bin")
nzx = read_write(file,"r",n=n);
nzx = reshape(nzx,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zy.bin")
nzy = read_write(file,"r",n=n);
nzy = reshape(nzy,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zz.bin")
nzz = read_write(file,"r",n=n);
nzz = reshape(nzz,n);

# ssa data
file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zx.bin")
szx = read_write(file,"r",n=n);
szx = reshape(szx,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zy.bin")
szy = read_write(file,"r",n=n);
szy = reshape(szy,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zz.bin")
szz = read_write(file,"r",n=n);
szz = reshape(szz,n);

# qssa data
file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zx.bin")
qzx = read_write(file,"r",n=n);
qzx = reshape(qzx,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zy.bin")
qzy = read_write(file,"r",n=n);
qzy = reshape(qzy,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zz.bin")
qzz = read_write(file,"r",n=n);
qzz = reshape(qzz,n);

# aqssa data
file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zx.bin")
azx = read_write(file,"r",n=n);
azx = reshape(azx,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zy.bin")
azy = read_write(file,"r",n=n);
azy = reshape(azy,n);

file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zz.bin")
azz = read_write(file,"r",n=n);
azz = reshape(azz,n);

function crg_mat_wigb(dx,dy,dz;
                      dt=0.004,
                      nt=100,n1=20,npanel=4,
                      irx=10,iry=10,
                      xtickpos = [16;42;67;93],
                      xticklab = ([x[5];x[10];x[15];x[20]] .- 1) .* 10,
                      ytickpos=1:50:99,
                      FontSize=16,
                      sc=0.5,
                      amx=1.0,
                      fname="noisy5d3c",
                      w=10,
                      h=15)
    t = collect((0:nt-1).*dt);

    tmpx = zeros(nt,n1,npanel);
    tmpy = zeros(nt,n1,npanel);
    tmpz = zeros(nt,n1,npanel);
    j = 1;
    for i in (5,10,15,20)
        tmpx[:,:,j] .= dx[:,:,i,10,10];
        tmpy[:,:,j] .= dy[:,:,i,10,10];
        tmpz[:,:,j] .= dz[:,:,i,10,10];
        j += 1;
    end

    tmp[1,:,:] .= 0;
    tt,ty,tx = size(tmp);
    x = zeros(tx*ty);

    c0=6;
    k=1; c=0;
    for ix = 1:tx
        c += c0;
        for iy = 1:ty
            x[k] = (ix-1)*(ty-1)+iy+c;
            k += 1;
        end
    end

    Dx = reshape(tmpx,nt,:);
    Dy = reshape(tmpy,nt,:);
    Dz = reshape(tmpz,nt,:);


    mat"close all"
    mat"h3 = tight_subplot(3,1,0.05,[0.1 0.07],[0.1 0.01])"
    mat"subplot(311)"
    mat"wigb($(Dx),$(sc),$x,$t,$(amx))"
    mat"title({'{(a)} Ux'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"ylabel({'t(s)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"set(gca,'TickLabelInterpreter','latex')"
    mat"set(gca,'FontSize',$FontSize)"
    mat"xticks($xtickpos)"
    mat"xticklabels($xticklab)"
    mat"yticks([0 0.2 0.396])"
    mat"yticklabels([0 0.2 0.4])"

    mat"subplot(312)"
    mat"wigb($(Dy),$(sc),$x,$t,$(amx))"
    mat"title({'{(b)} Uy'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"ylabel({'t(s)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"set(gca,'TickLabelInterpreter','latex')"
    mat"set(gca,'FontSize',$FontSize)"
    mat"xticks($xtickpos)"
    mat"xticklabels($xticklab)"
    mat"yticks([0 0.2 0.396])"
    mat"yticklabels([0 0.2 0.4])"

    mat"subplot(313)"
    mat"wigb($(Dz),$(sc),$x,$t,$(amx))"
    mat"title({'{(c)} Uz'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"ylabel({'t (s)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"xlabel({'sx (m)'},'FontSize',$FontSize,'Interpreter','latex')"
    mat"set(gca,'TickLabelInterpreter','latex')"
    mat"set(gca,'FontSize',$FontSize)"
    mat"xticks($xtickpos)"
    mat"xticklabels($xticklab)"
    mat"yticks([0 0.2 0.396])"
    mat"yticklabels([0 0.2 0.4])"

    mat"set(gcf,'PaperOrientation','landscape')"
    mat"set(gcf,'PaperSize',[$w $h])"
    mat"print(gcf,$(fname),'-dpdf','-fillpage')"
#    mat"exportgraphics(gcf,$(fname),'ContentType','vector')"
end
