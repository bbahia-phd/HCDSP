
# K = 1:2:50;
# kmax = length(K)
# o1  = zeros(eltype(dc),size(dc)...,kmax);
# o2  = zeros(eltype(dc),size(dc)...,kmax);
# o3  = zeros(eltype(dc),size(dc)...,kmax);
# tmp = zeros(eltype(dc),size(dc)...);

# rmax = 20;

# r1 = zeros(rmax,kmax);
# r2 = zeros(rmax,kmax);
# r3 = zeros(rmax,kmax);

# for k in 1:kmax
#     kk = K[k]

#     for r in 1:rmax
#         dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
#         INF = complex.(PadOp(dnx,nin=nin,npad=npad,flag="fwd"));
#         fft!(INF,1);

#         d  = copy(INF[iω,indx]);

#         tmp .= Op1(d,kk);
#         r1[r,k] = quality(tmp,dc);
#         #o1[:,:,kk] .= tmp;

#         tmp .= Op2(d,kk);
#         r2[r,k] = quality(tmp,dc);
#         #o2[:,:,kk] .= tmp;

#         tmp .= Op3(d,kk);
#         r3[r,k] = quality(tmp,dc);
#         #o3[:,:,kk] .= tmp;
#     end
#     @show [kk  mean(r1[:,k],dims=1) mean(r2[:,k],dims=1) mean(r3[:,k],dims=1)]
# end

# # Average
# r1r = vec(mean(r1,dims=1));
# r2r = vec(mean(r2,dims=1));
# r3r = vec(mean(r3,dims=1));

# std1 = vec(std(r1,dims=1));
# std2 = vec(std(r2,dims=1));
# std3 = vec(std(r3,dims=1));

# close("all");
# clf();

# errorbar(K,r1r,yerr=std1,fmt="-o",label="SVD")
# errorbar(K,r2r,yerr=std2,fmt="-o",label="rQR")
# errorbar(K,r3r,yerr=std3,fmt="-o",label="Lanczos")
# xlabel("'Rank' "*L"(k)",fontsize=15)
# ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
# legend()

# gcf()

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

# using HDF5

# fname = "./hcdsp_low_rank_compare_denoising"
# fid = h5open(fname, "cw")

# create_group(fid,"gains")
# fid["gains"]["svd"] = r1;
# fid["gains"]["rqr"] = r2;
# fid["gains"]["lanc"] = r3;

# create_group(fid,"gains/stats/mean")
# fid["gains/stats/mean"]["svd"]  = r1r
# fid["gains/stats/mean"]["rqr"]  = r2r
# fid["gains/stats/mean"]["lanc"] = r3r

# create_group(fid,"gains/stats/std")
# fid["gains/stats/std"]["svd"]  = std1
# fid["gains/stats/std"]["rqr"]  = std2
# fid["gains/stats/std"]["lanc"] = std3

# create_group(fid,"data")
# fid["data"]["clean"] = dzx
# fid["data"]["noisy"] = dnx

# close(fid)


cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using HCDSP
using PyPlot
using StatsBase,Statistics
using SeisMain, SeisPlot

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=100, ox2=0.0, dx2=10.0, nx2=100, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);
dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
dnx = decimate_traces(dnx,40);

nin = size(dzx); npad = nin;

INC = complex.(PadOp(dzx,nin=nin,npad=npad,flag="fwd"));
INF = complex.(PadOp(dnx,nin=nin,npad=npad,flag="fwd"));

fft!(INC,1);
fft!(INF,1);

fmin = 0.0; fmax = 50.0; dt = 0.004;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

iω = 15;

dc = copy(INC[iω,indx]);
dn = copy(INF[iω,indx]);
T = SamplingOp(dn);

imp_ssa(d,k) = HCDSP.imputation_op(d,T,SVDSSAOp,(k))
imp_rqr(d,k) = HCDSP.imputation_op(d,T,rQR,(k))
imp_ssa(d,k) = HCDSP.imputation_op(d,T,LANCSSAOp,(k))

dsvd = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,SVDSSAOp,(5))...);
dqr  = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,rQROp,(8))...);
dout = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,FSSAOp,(5))...);
