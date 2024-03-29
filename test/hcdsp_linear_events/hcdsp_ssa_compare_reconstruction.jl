cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using HCDSP
using PyPlot
using LinearAlgebra
using StatsBase,Statistics
using SeisMain, SeisPlot

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
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

# These act on a frequency slice d
imp_ssa(d,T,k) = HCDSP.imputation_op(d,T,SVDSSAOp,(k))
imp_rqr(d,T,k) = HCDSP.imputation_op(d,T,rQROp,(k))
imp_lan(d,T,k) = HCDSP.imputation_op(d,T,LANCSSAOp,(k))

# ranks to test
K = 1:2:50;
kmax = length(K);

# decimations to test
percs = 10:10:90;
pmax = length(percs);

# for outputs
tmp = zeros(eltype(dc),size(dc)...);
# o1  = zeros(eltype(dc),size(dc)...,pmax,kmax);
# o2  = zeros(eltype(dc),size(dc)...,pmax,kmax);
# o3  = zeros(eltype(dc),size(dc)...,pmax,kmax);

rmax = 20;

r1 = zeros(rmax,pmax,kmax);
r2 = zeros(rmax,pmax,kmax);
r3 = zeros(rmax,pmax,kmax);

for k in 1:kmax
    kk = K[k]

    for p in 1:pmax
        perc = percs[p]

        for r in 1:rmax
            dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
            dnx .= decimate_traces(dnx,perc);
            
            INF = complex.(PadOp(dnx,nin=nin,npad=npad,flag="fwd"));
            fft!(INF,1);
            d  = copy(INF[iω,indx]);

            T = SamplingOp(d);

            tmp .= imp_ssa(d,T,kk);
            r1[r,p,k] = quality(tmp,dc);
            #o1[:,:,kk] .= tmp;
    
            tmp .= imp_rqr(d,T,kk);
            r2[r,p,k] = quality(tmp,dc);
            #o2[:,:,kk] .= tmp;
    
            tmp .= imp_lan(d,T,kk);
            r3[r,p,k] = quality(tmp,dc);
            #o3[:,:,kk] .= tmp;
        end       
        @show [kk perc mean(r1[:,p,k],dims=1) mean(r2[:,p,k],dims=1) mean(r3[:,p,k],dims=1)]    
    end
end

using HDF5

fname = joinpath(homedir(),"projects/HCDSP/hcdsp_ssa_reconstruction.h5")
fid = h5open(fname, "w")

create_group(fid,"gains")
fid["gains"]["svd"] = r1;
fid["gains"]["rqr"] = r2;
fid["gains"]["lanc"] = r3;

close(fid)

# Calls for different SSA-based reconstruction
# dsvd = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,SVDSSAOp,(5))...);
# dqr  = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,rQROp,(20))...);
# dqr  = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,QRFSSAOp,(5))...);
# dout = fx_process(dnx,dt,fmin,fmax,HCDSP.imputation_op,(T,FSSAOp,(5))...);

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

