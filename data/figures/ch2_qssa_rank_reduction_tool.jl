cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using HDF5
using HCDSP
using PyPlot
using LinearAlgebra
using StatsBase,Statistics
using SeisMain, SeisPlot

fname = joinpath(homedir(),"projects/HCDSP/data/hcdsp_qssa_compare_reconstruction.h5")

r1 = h5read(fname, "gains/svd")
r2 = h5read(fname, "gains/rqr")
r3 = h5read(fname, "gains/lanc")

# Average
r1r = mean(r1,dims=1);
r2r = mean(r2,dims=1);
r3r = mean(r3,dims=1);

std1 = std(r1,dims=1);
std2 = std(r2,dims=1);
std3 = std(r3,dims=1);

# ranks to test
K = 1:2:50;
kmax = length(K);

# decimations to test
percs = 10:10:90;
pmax = length(percs);
p = 6; perc = percs[p];

close("all");clf();

errorbar(K,r1r[1,p,:],yerr=std1[1,p,:],fmt="-o",label="SVD")
errorbar(K,r2r[1,p,:],yerr=std2[1,p,:],fmt="-o",label="rQR")
errorbar(K,r3r[1,p,:],yerr=std3[1,p,:],fmt="-o",label="Lanczos")
xlabel("'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
legend()

gcf()

# close("all");
# clf();
# plot(K,r1r);
# plot(K,r2r);
# plot(K,r3r);
# gcf()

# close("all");clf()
# fignum="time_domain_data";
# figure(fignum,figsize=(9,5));

# subplot(1,2,1);
# SeisPlotTX(dzx[:,1:2:80,1],
#         dy=dt,
#         dx=10,
#         cmap="gray",
#         style="overlay",
#         fignum=fignum,
#         xcur=2,
#         title="(a)",
#         ylabel="Time (s)",
#         xlabel="Offset (m)")

# subplot(1,2,2);
# SeisPlotTX(dnx[:,1:2:80,1],
#         dy=dt,
#         dx=10,
#         cmap="gray",
#         style="overlay",
#         fignum=fignum,
#         xcur=2,
#         title="(b)",
#         ylabel="Time (s)",
#         xlabel="Offset (m)")
# tight_layout();
# gcf()

# close(fid)