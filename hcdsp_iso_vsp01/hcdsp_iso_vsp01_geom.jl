cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate("./dev/hcdsp/.")
Pkg.status()

using Revise

using FFTW
using PyPlot
using SeisMain, SeisPlot
using HCDSP,SNR,IterativeMethods

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz_crg1350.seis"));

dt = ext_zz.d1;
nt = ext_zz.n1;

ntr = ext_zz.n2;

rz_init = 1350.0f0; # m
rz_end  = 1850.01f0; # m
drz = 16.667 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr = length(rz_axis)

nsline = 205; # number of source lines
ns = 205; # number of sources within lines

@assert ns*nsline == ntr == size(dzz,2)

dzx = reshape(dzx,(nt,ns,nsline));
dzy = reshape(dzy,(nt,ns,nsline));
dzz = reshape(dzz,(nt,ns,nsline));

j = 31;
comp = [dzx[:,:,j] dzy[:,:,j] dzz[:,:,j]];

SeisPlotTX(comp,
           cmap="gray",
           pclip=90,
           wbox=9,
           hbox=9,
           oy=ext_zx.o1,
           dy=ext_zx.d1);
gcf()
j = 50;
comp = [dzx[:,j,:] dzy[:,j,:] dzz[:,j,:]];

SeisPlotTX(comp,
           cmap="gray",
           pclip=90,
           wbox=9,
           hbox=9,
           oy=ext_zx.o1,
           dy=ext_zx.d1);
gcf()

j = 120;
comp = [dzx[j,:,:] dzy[j,:,:] dzz[j,:,:]];

SeisPlotTX(comp,
           cmap="gray",
           pclip=90,
           wbox=9,
           hbox=4);
gcf()


Q = quaternion(dzx,dzy,dzz);

psize = (128,32,32);
polap = (50,50,50);
smin  = (1,1,1);
smax  = size(Q)

qpatch,pid = fwdPatchOp(Q,psize,polap,smin,smax);

ax = quaternion(1.0f0,0.0f0,0.0f0)

side = "left"

Qf = qfft(qpatch[15][:,:,10],ax,side,1);

clf();
plot(imagi.(Qf[30,:]));gcf()