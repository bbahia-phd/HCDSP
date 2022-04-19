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

# Component-wise denoising
drx = fx_process(dnx, 0.004, 0, 100, SVDSSAOp, (2));

# prediction_quality
Rx = quality(drx,dzx)