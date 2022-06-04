pwd()

using Pkg
Pkg.activate("./../")
Pkg.status()

using Revise
using HDF5
using PyPlot
using StatsBase,Statistics


fname="hcdsp_reconstruct_mixed_linear_events.h5"
fid=h5open(fname,"r")

# extract gains
ax = read(fid["gains"]["aquater"]["x"])[:,:,1:11];
ay = read(fid["gains"]["aquater"]["y"])[:,:,1:11];
az = read(fid["gains"]["aquater"]["z"])[:,:,1:11];

qx = read(fid["gains"]["quater"]["x"])[:,:,1:11];
qy = read(fid["gains"]["quater"]["y"])[:,:,1:11];
qz = read(fid["gains"]["quater"]["z"])[:,:,1:11];

sx = read(fid["gains"]["real"]["x"])[:,:,1:11];
sy = read(fid["gains"]["real"]["y"])[:,:,1:11];
sz = read(fid["gains"]["real"]["z"])[:,:,1:11];

close(fid)

fix_mean(x) = dropdims(mean(x,dims=1),dims=1)
fix_std(x)  = dropdims(std(x,dims=1),dims=1)

ax_m = fix_mean(ax);#dropdims(mean(qx,dims=1),dims=1)
ay_m = fix_mean(ay);#dropdims(mean(qy,dims=1),dims=1)
az_m = fix_mean(az);#dropdims(mean(qz,dims=1),dims=1)

qx_m = fix_mean(qx);#dropdims(mean(qx,dims=1),dims=1)
qy_m = fix_mean(qy);#dropdims(mean(qy,dims=1),dims=1)
qz_m = fix_mean(qz);#dropdims(mean(qz,dims=1),dims=1)

sx_m = fix_mean(sx);#dropdims(mean(sx,dims=1),dims=1)
sy_m = fix_mean(sy);#dropdims(mean(sy,dims=1),dims=1)
sz_m = fix_mean(sz);#dropdims(mean(sz,dims=1),dims=1)

ax_s = fix_std(ax);#dropdims(std(qx,dims=1),dims=1)
ay_s = fix_std(ay);#dropdims(std(qy,dims=1),dims=1)
az_s = fix_std(az);#dropdims(std(qz,dims=1),dims=1)

qx_s = fix_std(qx);#dropdims(std(qx,dims=1),dims=1)
qy_s = fix_std(qy);#dropdims(std(qy,dims=1),dims=1)
qz_s = fix_std(qz);#dropdims(std(qz,dims=1),dims=1)

sx_s = fix_std(sx);#dropdims(std(sx,dims=1),dims=1)
sy_s = fix_std(sy);#dropdims(std(sy,dims=1),dims=1)
sz_s = fix_std(sz);#dropdims(std(sz,dims=1),dims=1)

K=1:2:21;

perc=10:10:90;

p=6;
j=2;

fig_name="rank_red_compare_ssa_perc60"

lsvd="SSA"
lrqr="QSSA"
llan="AQSSA"

close("all");
figure(fig_name,figsize=(12,4))

subplot(1,3,1)
errorbar(perc,sx_m[:,p],yerr=sx_s[:,p],fmt="-o",label=lsvd,color="blue")
errorbar(perc,qx_m[:,p],yerr=qx_s[:,p],fmt="-o",label=lrqr,color="red")
errorbar(perc,ax_m[:,p+j],yerr=ax_s[:,p+j],fmt="-o",label=llan,color="black")
ylim(10,60)
xlabel(L"Decimation "*L"(\%)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(a) "*L" R_{x}")
legend()

subplot(1,3,2)
errorbar(perc,sy_m[:,p],yerr=sy_s[:,p],fmt="-o",label=lsvd,color="blue")
errorbar(perc,qy_m[:,p],yerr=qy_s[:,p],fmt="-o",label=lrqr,color="red")
errorbar(perc,ay_m[:,p+j],yerr=ay_s[:,p+j],fmt="-o",label=llan,color="black")
ylim(10,60)
xlabel(L"Decimation "*L"(\%)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(b) "*L" R_{y}")
legend()

subplot(1,3,3)
errorbar(perc,sz_m[:,p],yerr=sz_s[:,p],fmt="-o",label=lsvd,color="blue")
errorbar(perc,qz_m[:,p],yerr=qz_s[:,p],fmt="-o",label=lrqr,color="red")
errorbar(perc,az_m[:,p+j],yerr=az_s[:,p+j],fmt="-o",label=llan,color="black")
ylim(10,60)
xlabel(L"Decimation "*L"(\%)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(c) "*L" R_{z}")

legend()
tight_layout() 

savefig(fig_name)


fig_name="rank_red_compare_ssa_rank60"

lsvd="SSA"
lrqr="QSSA"
llan="AQSSA"

close("all");
figure(fig_name,figsize=(12,4))

subplot(1,3,1)
errorbar(K,sx_m[p,:],yerr=sx_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qx_m[p,:],yerr=qx_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,ax_m[p,:],yerr=ax_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,55)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(a) "*L" R_{x}")
legend()

subplot(1,3,2)
errorbar(K,sy_m[p,:],yerr=sy_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qy_m[p,:],yerr=qy_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,ay_m[p,:],yerr=ay_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,55)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(b) "*L" R_{y}")
legend()

subplot(1,3,3)
errorbar(K,sz_m[p,:],yerr=sz_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qz_m[p,:],yerr=qz_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,az_m[p,:],yerr=az_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,55)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(c) "*L" R_{z}")

legend()
tight_layout() 

savefig(fig_name)
