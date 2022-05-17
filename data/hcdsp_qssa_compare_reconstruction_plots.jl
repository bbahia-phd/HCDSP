pwd()

using Pkg
Pkg.activate("./../")
Pkg.status()

using Revise
using HDF5
using PyPlot
using StatsBase,Statistics


fname="./hcdsp_qssa_compare_reconstruction.h5"
fid=h5open(fname,"r")

# extract gains
lx = read(fid["gains"]["lanc"]["x"]);
ly = read(fid["gains"]["lanc"]["y"]);
lz = read(fid["gains"]["lanc"]["z"]);

qx = read(fid["gains"]["rqr"]["x"]);
qy = read(fid["gains"]["rqr"]["y"]);
qz = read(fid["gains"]["rqr"]["z"]);

sx = read(fid["gains"]["svd"]["x"]);
sy = read(fid["gains"]["svd"]["y"]);
sz = read(fid["gains"]["svd"]["z"]);

close(fid)

fix_mean(x) = reverse(dropdims(mean(x,dims=1),dims=1),dims=1)
fix_std(x) =  reverse(dropdims(std(x,dims=1),dims=1),dims=1)

lx_m = fix_mean(lx);#dropdims(mean(lx,dims=1),dims=1)
ly_m = fix_mean(ly);#dropdims(mean(ly,dims=1),dims=1)
lz_m = fix_mean(lz);#dropdims(mean(lz,dims=1),dims=1)

qx_m = fix_mean(qx);#dropdims(mean(qx,dims=1),dims=1)
qy_m = fix_mean(qy);#dropdims(mean(qy,dims=1),dims=1)
qz_m = fix_mean(qz);#dropdims(mean(qz,dims=1),dims=1)

sx_m = fix_mean(sx);#dropdims(mean(sx,dims=1),dims=1)
sy_m = fix_mean(sy);#dropdims(mean(sy,dims=1),dims=1)
sz_m = fix_mean(sz);#dropdims(mean(sz,dims=1),dims=1)

lx_s = fix_std(lx);#dropdims(std(lx,dims=1),dims=1)
ly_s = fix_std(ly);#dropdims(std(ly,dims=1),dims=1)
lz_s = fix_std(lz);#dropdims(std(lz,dims=1),dims=1)

qx_s = fix_std(qx);#dropdims(std(qx,dims=1),dims=1)
qy_s = fix_std(qy);#dropdims(std(qy,dims=1),dims=1)
qz_s = fix_std(qz);#dropdims(std(qz,dims=1),dims=1)

sx_s = fix_std(sx);#dropdims(std(sx,dims=1),dims=1)
sy_s = fix_std(sy);#dropdims(std(sy,dims=1),dims=1)
sz_s = fix_std(sz);#dropdims(std(sz,dims=1),dims=1)

K=1:2:50;

perc=10:10:90;

p=5;

fig_name="rank_red_tool_error_bar"

lsvd="SVD"
lrqr="rQR"
llan="Lanczos"

close("all");
figure(fig_name,figsize=(12,4))

subplot(1,3,1)
errorbar(K,sx_m[p,:],yerr=sx_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qx_m[p,:],yerr=lx_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,lx_m[p,:],yerr=qx_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(a) "*L" R_{x}")
legend()

subplot(1,3,2)
errorbar(K,sy_m[p,:],yerr=sy_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qy_m[p,:],yerr=ly_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,ly_m[p,:],yerr=qy_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(b) "*L" R_{y}")
legend()

subplot(1,3,3)
errorbar(K,sz_m[p,:],yerr=sz_s[p,:],fmt="-o",label=lsvd,color="blue")
errorbar(K,qz_m[p,:],yerr=lz_s[p,:],fmt="-o",label=lrqr,color="red")
errorbar(K,lz_m[p,:],yerr=qz_s[p,:],fmt="-o",label=llan,color="black")
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
title("(c) "*L" R_{z}")

legend()
tight_layout() 

fig_name="rank_red_tool_error_bar"

lsvd="SVD"
lrqr="rQR"
llan="Lanczos"

perc=string.(reverse(10:10:90))[1:2:end];
k = string.(K[1:2:end]);
 
close("all");
figure(fig_name,figsize=(12,12))

subplot(331)
imshow(sx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(a) SVD "*L" R_{x}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(332)
imshow(sy_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(b) SVD "*L" R_{y}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(333)
imshow(sz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(c) SVD "*L" R_{z}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(334)
imshow(lx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(d) Lanczos "*L" R_{x}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(335)
imshow(ly_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(e) Lanczos "*L" R_{y}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(336)
imshow(lz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(f) Lanczos "*L" R_{z}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(337)
imshow(qx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(g) rQR "*L" R_{x}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(338)
imshow(qy_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(h) rQR "*L" R_{y}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

subplot(339)
imshow(qz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")
title("(i) rQR "*L" R_{z}")
xlabel(L"k")
ylabel("Decimation (%)")
ax=gca()
ax.set_xticks(0:2:24)
ax.set_xticklabels(k)
ax.set_yticks(0:2:8)
ax.set_yticklabels(perc)
colorbar()

tight_layout()

#=
fig_name="rank_red_tool"

figure(fig_name,fig_size=(9,3))

subplot(1,3,1)
plot(K,sx_m[p,:])
plot(K,lx_m[p,:])
plot(K,qx_m[p,:])

subplot(1,3,2)
plot(K,sy_m[p,:])
plot(K,ly_m[p,:])
plot(K,qy_m[p,:])

subplot(1,3,3)
plot(K,sz_m[p,:])
plot(K,lz_m[p,:])
plot(K,q_m[p,:])

xlabel("'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
legend()
=#
