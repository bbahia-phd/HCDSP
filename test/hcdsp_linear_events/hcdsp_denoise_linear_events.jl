cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using PyPlot
using SeisMain, SeisPlot
using HCDSP,IterativeMethods

params_zx = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);
dnx = SeisAddNoise(dzx, 0.5, db=true, L=9);

params_zy = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[-1.0,1.0], f0=20.0)
dzy = SeisLinearEvents(; params_zy...);
dny = SeisAddNoise(dzy, 1.0, db=true, L=9);

params_zz = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[-1.0,-1.0], f0=20.0)
dzz = SeisLinearEvents(; params_zz...);
dnz = SeisAddNoise(dzz, 3.0, db=true, L=9);

# Define time domain quaternion
Q = quaternion(dnx,dny,dnz);

# Quaternion denoising
QQ = fx_process(Q, 0.004, 0, 100, SVDSSAOp, (4));

# Component-wise denoising
drx = fx_process(dnx, 0.004, 0, 100, SVDSSAOp, (2));
dry = fx_process(dny, 0.004, 0, 100, SVDSSAOp, (2));
drz = fx_process(dnz, 0.004, 0, 100, SVDSSAOp, (2));

# prediction_quality

Rx = quality(drx,dzx)
Ry = quality(dry,dzy)
Rz = quality(drz,dzz)

Qx = quality(imagi.(QQ),dzx)
Qy = quality(imagj.(QQ),dzy)
Qz = quality(imagk.(QQ),dzz)

dqx = imagi.(QQ);
dqy = imagj.(QQ);
dqz = imagk.(QQ);

# get slice j
j = 20;

figure("panel",figsize=(10,10));

subplot(3,1,1);
SeisPlotTX([dzx[:,:,j] dnx[:,:,j] drx[:,:,j] dqx[:,:,j]],fignum="panel",style="overlay");

subplot(3,1,2);
SeisPlotTX([dzy[:,:,j] dny[:,:,j] dry[:,:,j] dqy[:,:,j]],fignum="panel",style="overlay");

subplot(3,1,3);
SeisPlotTX([dzz[:,:,j] dnz[:,:,j] drz[:,:,j] dqz[:,:,j]],fignum="panel",style="overlay");

gcf()