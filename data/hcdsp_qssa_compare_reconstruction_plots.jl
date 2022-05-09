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

lx_m = dropdims(mean(lx,dims=1),dims=1)
ly_m = dropdims(mean(ly,dims=1),dims=1)
lz_m = dropdims(mean(lz,dims=1),dims=1)

qx_m = dropdims(mean(qx,dims=1),dims=1)
qy_m = dropdims(mean(qy,dims=1),dims=1)
qz_m = dropdims(mean(qz,dims=1),dims=1)

sx_m = dropdims(mean(sx,dims=1),dims=1)
sy_m = dropdims(mean(sy,dims=1),dims=1)
sz_m = dropdims(mean(sz,dims=1),dims=1)

lx_s = dropdims(std(lx,dims=1),dims=1)
ly_s = dropdims(std(ly,dims=1),dims=1)
lz_s = dropdims(std(lz,dims=1),dims=1)

qx_s = dropdims(std(qx,dims=1),dims=1)
qy_s = dropdims(std(qy,dims=1),dims=1)
qz_s = dropdims(std(qz,dims=1),dims=1)

sx_s = dropdims(std(sx,dims=1),dims=1)
sy_s = dropdims(std(sy,dims=1),dims=1)
sz_s = dropdims(std(sz,dims=1),dims=1)

K=1:2:50;

perc=10:10:90;

p=5;

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
plot(K,qz_m[p,:])

fig_name="rank_red_tool_error_bar"

lsvd="SVD"
lrqr="rQR"
llan="Lanczos"

close("all");
figure(fig_name,figsize=(9,3))

subplot(1,3,1)
errorbar(K,sx_m[p,:],yerr=sx_s[p,:],fmt="-o",label=lsvd)
errorbar(K,lx_m[p,:],yerr=lx_s[p,:],fmt="-o",label=lrqr)
errorbar(K,qx_m[p,:],yerr=qx_s[p,:],fmt="-o",label=llan)
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)

subplot(1,3,2)
errorbar(K,sy_m[p,:],yerr=sy_s[p,:],fmt="-o",label=lsvd)
errorbar(K,ly_m[p,:],yerr=ly_s[p,:],fmt="-o",label=lrqr)
errorbar(K,qy_m[p,:],yerr=qy_s[p,:],fmt="-o",label=llan)
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)

subplot(1,3,3)
errorbar(K,sz_m[p,:],yerr=sz_s[p,:],fmt="-o",label=lsvd)
errorbar(K,lz_m[p,:],yerr=lz_s[p,:],fmt="-o",label=lrqr)
errorbar(K,qz_m[p,:],yerr=qz_s[p,:],fmt="-o",label=llan)
ylim(0,60)
xlabel(L"'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)

tight_layout() 


subplot(331)
imshow(sx_m,aspect="auto",vmin=0,vmax=50)

subplot(332)
imshow(lx_m,aspect="auto",vmin=0,vmax=50)

subplot(333)
imshow(qx_m,aspect="auto",vmin=0,vmax=50)

subplot(334)
imshow(sy_m,aspect="auto",vmin=0,vmax=50)

subplot(335)
imshow(ly_m,aspect="auto",vmin=0,vmax=50)

subplot(336)
imshow(qy_m,aspect="auto",vmin=0,vmax=50)

subplot(337)
imshow(sz_m,aspect="auto",vmin=0,vmax=50)

subplot(338)
imshow(lz_m,aspect="auto",vmin=0,vmax=50)

subplot(339)
imshow(qz_m,aspect="auto",vmin=0,vmax=50)



tight_layout()


subplot(331)
imshow(sx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian",cmap="seismic")

subplot(332)
imshow(lx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(333)
imshow(qx_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(334)
imshow(sy_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(335)
imshow(ly_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(336)
imshow(qy_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(337)
imshow(sz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(338)
imshow(lz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

subplot(339)
imshow(qz_m,aspect="auto",vmin=0,vmax=50,interpolation="gaussian")

tight_layout()
