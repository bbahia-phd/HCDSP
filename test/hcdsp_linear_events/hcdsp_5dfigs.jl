# pwd()

# dev_dir = "/home/bbahia/projects";
# dev_dir = "/home/brenobahia/projects";
# cd(dev_dir)
# pwd()

# using Pkg
# Pkg.activate(joinpath(dev_dir,"HCDSP"))
# Pkg.status()

# using Revise
# using LinearAlgebra
# using FFTW
# using HCDSP

# using PyPlot
# using SeisMain, SeisPlot

# using MAT,MATLAB
# mat"addpath('/home/brenobahia/projects/QSEIS/src/Plotting/')"
# #mat"addpath('/home/bbahia/projects/QSEIS/src/Plotting/')"

# using HDF5

# # data dir home
# #dir_path  = "/media/bbahia/DATA/seismic_data/linear5d";
# dir_path  = joinpath(dev_dir,"files/linear5d");

# n  = nt,n1,n2,n3,n4 = 100,20,20,20,20;
# dt = 0.004;

# # ideal data
# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zx.bin")
# dzx = read_write(file,"r",n=n);
# dzx = reshape(dzx,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zy.bin")
# dzy = read_write(file,"r",n=n);
# dzy = reshape(dzy,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_zz.bin")
# dzz = read_write(file,"r",n=n);
# dzz = reshape(dzz,n);

# # noisy data
# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zx.bin")
# nzx = read_write(file,"r",n=n);
# nzx = reshape(nzx,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zy.bin")
# nzy = read_write(file,"r",n=n);
# nzy = reshape(nzy,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_noisy_zz.bin")
# nzz = read_write(file,"r",n=n);
# nzz = reshape(nzz,n);

# # ssa data
# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zx.bin")
# szx = read_write(file,"r",n=n);
# szx = reshape(szx,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zy.bin")
# szy = read_write(file,"r",n=n);
# szy = reshape(szy,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_ssa_zz.bin")
# szz = read_write(file,"r",n=n);
# szz = reshape(szz,n);

# # qssa data
# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zx.bin")
# qzx = read_write(file,"r",n=n);
# qzx = reshape(qzx,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zy.bin")
# qzy = read_write(file,"r",n=n);
# qzy = reshape(qzy,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_qssa_zz.bin")
# qzz = read_write(file,"r",n=n);
# qzz = reshape(qzz,n);

# # aqssa data
# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zx.bin")
# azx = read_write(file,"r",n=n);
# azx = reshape(azx,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zy.bin")
# azy = read_write(file,"r",n=n);
# azy = reshape(azy,n);

# file = joinpath(dir_path,"hcdsp_recon_3d3c_linear_events_aqssa_zz.bin")
# azz = read_write(file,"r",n=n);
# azz = reshape(azz,n);

close("all");

w=6; h=10;

indx_to_plot = 1:2:n1;

fname="my_fig";
plot_param = Dict(:fignum => fname,
                  :style  => "overlay",
                  :xcur   => 2.5,
                  :cmap   => "gray",
                  :xticks => [5] );

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(dzx[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax12)
SeisPlotTX(dzx[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(dzx[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(dzy[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax22)
SeisPlotTX(dzy[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(dzy[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(dzz[:,indx_to_plot,5,1,1];plot_param...,xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax32)
SeisPlotTX(dzz[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z");
subplot(ax33)
SeisPlotTX(dzz[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

tight_layout()
gcf()

savefig("./clean_input")

close("all");

w=6; h=10;

indx_to_plot = 1:2:n1;

fname="my_fig";
plot_param = Dict(:fignum => fname,
                  :style  => "overlay",
                  :xcur   => 2.5,
                  :cmap   => "gray",
                  :xticks => [5] );

figure(fname,figsize=(w,h))
fig = gcf()

ax11 = subplot2grid((3,3),(0,0),1,1)
ax12 = subplot2grid((3,3),(0,1),1,1)
ax13 = subplot2grid((3,3),(0,2),1,1)

subplot(ax11)
SeisPlotTX(nzx[:,indx_to_plot,5,1,1];plot_param..., xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax12)
SeisPlotTX(nzx[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_x");
subplot(ax13)
SeisPlotTX(nzx[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

ax21 = subplot2grid((3,3),(1,0),1,1)
ax22 = subplot2grid((3,3),(1,1),1,1)
ax23 = subplot2grid((3,3),(1,2),1,1)

subplot(ax21)
SeisPlotTX(nzy[:,indx_to_plot,5,1,1];plot_param..., xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax22)
SeisPlotTX(nzy[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_y");
subplot(ax23)
SeisPlotTX(nzy[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

ax31 = subplot2grid((3,3),(2,0),1,1)
ax32 = subplot2grid((3,3),(2,1),1,1)
ax33 = subplot2grid((3,3),(2,2),1,1)

subplot(ax31)
SeisPlotTX(nzz[:,indx_to_plot,5,1,1];plot_param..., xticklabels=[L"$s_x = 50$ (m) "]);
subplot(ax32)
SeisPlotTX(nzz[:,indx_to_plot,10,1,1];plot_param...,xticklabels=[L"$s_x = 100$ (m) "],title=L"{\bf U}_z");
subplot(ax33)
SeisPlotTX(nzz[:,indx_to_plot,20,1,1];plot_param...,xticklabels=[L"$s_x = 200$ (m) "]);

tight_layout()
gcf()

savefig("./noisy_input")
