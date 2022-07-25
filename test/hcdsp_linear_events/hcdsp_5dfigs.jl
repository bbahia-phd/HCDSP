 pwd()

 dev_dir = "/home/bbahia/projects";
# dev_dir = "/home/brenobahia/projects";
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
 #mat"addpath('/home/brenobahia/projects/QSEIS/src/Plotting/')"
 mat"addpath('/home/bbahia/projects/QSEIS/src/Plotting/')"

 using HDF5

# data dir home
dir_path  = "/media/bbahia/DATA/seismic_data/linear5d";
# dir_path  = joinpath(dev_dir,"files/linear5d");

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

close("all");

w=4; h=10;

indx_to_plot = 1:2:n1;

fname="my_fig";
plot_param = Dict(:fignum => fname,
                  :style  => "overlay",
                  :xcur   => 2.5,
                  :cmap   => "gray",
                  :xticks => [],
                  :scal   => 0.5,
                  :oy     => 0.0,
                  :dy     => dt,
                  :ylabel => "Time",
                  :yunits => "(s)",
                  :yticklabels => ["0.0","0.2","0.4"],
                  :yticks      => [0,0.2,0.396],
                  :labelsize   => 10);

plot_param_1st = Dict(:fignum => fname,
                      :style  => "overlay",
                      :xcur   => 2.5,
                      :cmap   => "gray",
                      :xticks => [],
                      :scal   => 0.5,
                      :oy     => 0.0,
                      :dy     => dt,
                      :ylabel => "",
                      :yunits => "",
                      :yticklabels => [],
                      :yticks      => [],
                      :labelsize   => 10);

plot_param_2nd = Dict(:fignum => fname,
                      :style  => "overlay",
                      :xcur   => 2.5,
                      :cmap   => "gray",
                      :scal   => 0.5,
                      :oy     => 0.0,
                      :dy     => dt,
                      :ylabel => "",
                      :yunits => "",
                      :yticklabels => [],
                      :yticks      => [],
                      :labelsize   => 10);

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(dzx[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax12)          
SeisPlotTX(dzx[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(dzx[:,indx_to_plot,20,1,1];plot_param_1st...);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(dzy[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax22)
SeisPlotTX(dzy[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(dzy[:,indx_to_plot,20,1,1];plot_param_1st...);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(dzz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "],xticks=[5]);
subplot(ax32)
SeisPlotTX(dzz[:,indx_to_plot,10,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z",xticks=[5]);
subplot(ax33)
SeisPlotTX(dzz[:,indx_to_plot,20,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 200$ (m) "],xticks=[5]);

subplots_adjust(wspace=0.05,hspace=0.15,left=0.2)
gcf()

savefig("./clean_input")

#tight_layout()
#SeisPlotTX(dzx[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "]);
#SeisPlotTX(dzx[:,indx_to_plot,10,1,1];plot_param_1st...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_x");
#SeisPlotTX(dzx[:,indx_to_plot,20,1,1];plot_param_1st...,xticklabels=[L"$s_x = 200$ (m) "]);

close("all");

w=4; h=10;

indx_to_plot = 1:2:n1;

fname="my_fig";
plot_param = Dict(:fignum => fname,
                  :style  => "overlay",
                  :xcur   => 2.5,
                  :cmap   => "gray",
                  :xticks => [],
                  :scal   => 0.5,
                  :oy     => 0.0,
                  :dy     => dt,
                  :ylabel => "Time",
                  :yunits => "(s)",
                  :yticklabels => ["0.0","0.2","0.4"],
                  :yticks      => [0,0.2,0.396],
                  :labelsize   => 10);

plot_param_1st = Dict(:fignum => fname,
                      :style  => "overlay",
                      :xcur   => 2.5,
                      :cmap   => "gray",
                      :xticks => [],
                      :scal   => 0.5,
                      :oy     => 0.0,
                      :dy     => dt,
                      :ylabel => "",
                      :yunits => "",
                      :yticklabels => [],
                      :yticks      => [],
                      :labelsize   => 10);

plot_param_2nd = Dict(:fignum => fname,
                      :style  => "overlay",
                      :xcur   => 2.5,
                      :cmap   => "gray",
                      :scal   => 0.5,
                      :oy     => 0.0,
                      :dy     => dt,
                      :ylabel => "",
                      :yunits => "",
                      :yticklabels => [],
                      :yticks      => [],
                      :labelsize   => 10);

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(nzx[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax12)          
SeisPlotTX(nzx[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(nzx[:,indx_to_plot,20,1,1];plot_param_1st...);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(nzy[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax22)
SeisPlotTX(nzy[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(nzy[:,indx_to_plot,20,1,1];plot_param_1st...);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(nzz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "],xticks=[5]);
subplot(ax32)
SeisPlotTX(nzz[:,indx_to_plot,10,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z",xticks=[5]);
subplot(ax33)
SeisPlotTX(nzz[:,indx_to_plot,20,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 200$ (m) "],xticks=[5]);

subplots_adjust(wspace=0.05,hspace=0.15,left=0.2)
gcf()

savefig("./noisy_input")

<<<<<<< Updated upstream
dev_dir = "/dev/Breno_GOM/projects";
#dev_dir = "/home/bbahia/projects";
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
dir_path = joinpath(dev_dir,"files/linear_5d");
#dir_path  = "/media/bbahia/DATA/seismic_data/linear5d";

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

    tmpx[1,:,:] .= 0;
    tmpy[1,:,:] .= 0;
    tmpz[1,:,:] .= 0;
    tt,ty,tx = size(tmpx);
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
=======
close("all");

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(szx[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax12)          
SeisPlotTX(szx[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(szx[:,indx_to_plot,20,1,1];plot_param_1st...);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(szy[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax22)
SeisPlotTX(szy[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(szy[:,indx_to_plot,20,1,1];plot_param_1st...);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(szz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "],xticks=[5]);
subplot(ax32)
SeisPlotTX(szz[:,indx_to_plot,10,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z",xticks=[5]);
subplot(ax33)
SeisPlotTX(szz[:,indx_to_plot,20,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 200$ (m) "],xticks=[5]);

subplots_adjust(wspace=0.05,hspace=0.15,left=0.2)
gcf()

savefig("./ssa_output")

close("all");

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(qzx[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax12)          
SeisPlotTX(qzx[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(qzx[:,indx_to_plot,20,1,1];plot_param_1st...);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(qzy[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax22)
SeisPlotTX(qzy[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(qzy[:,indx_to_plot,20,1,1];plot_param_1st...);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(qzz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "],xticks=[5]);
subplot(ax32)
SeisPlotTX(qzz[:,indx_to_plot,10,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z",xticks=[5]);
subplot(ax33)
SeisPlotTX(qzz[:,indx_to_plot,20,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 200$ (m) "],xticks=[5]);

subplots_adjust(wspace=0.05,hspace=0.15,left=0.2)
gcf()

savefig("./qssa_output") 

close("all");

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(azx[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax12)          
SeisPlotTX(azx[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(azx[:,indx_to_plot,20,1,1];plot_param_1st...);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(azy[:,indx_to_plot,5,1,1];plot_param...);
subplot(ax22)
SeisPlotTX(azy[:,indx_to_plot,10,1,1];plot_param_1st...,title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(azy[:,indx_to_plot,20,1,1];plot_param_1st...);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(azz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "],xticks=[5]);
subplot(ax32)
SeisPlotTX(azz[:,indx_to_plot,10,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z",xticks=[5]);
subplot(ax33)
SeisPlotTX(azz[:,indx_to_plot,20,1,1];plot_param_2nd...,xticklabels=[L"$s_x = 200$ (m) "],xticks=[5]);

subplots_adjust(wspace=0.05,hspace=0.15,left=0.2)
gcf()

savefig("./aqssa_output") 
>>>>>>> Stashed changes
