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
nx1=100, ox2=0.0, dx2=10.0, nx2=1, ox3=0.0, dx3=10.0,
nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);

params_zy = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
nx1=100, ox2=0.0, dx2=10.0, nx2=1, ox3=0.0, dx3=10.0,
nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
amp=[-1.0,1.0], f0=20.0)
dzy = SeisLinearEvents(; params_zy...);

params_zz = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
nx1=100, ox2=0.0, dx2=10.0, nx2=1, ox3=0.0, dx3=10.0,
nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
amp=[-1.0,-1.0], f0=20.0)
dzz = SeisLinearEvents(; params_zz...);

Q = quaternion(dzx,dzy,dzz);

ax = quaternion(1.0,0.0,0.0)

side = "left"

Qf = qfft(Q,ax,side,1);

qs = Qf[50,:];

H = HankelOp(qs) # Adjoint operator A'(S'(X))

U,S,V = tsvd(H,8);

