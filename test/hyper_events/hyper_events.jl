pwd()
cd(joinpath(homedir(),"projects/HCDSP/test/hyper_events"))

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using SeisMain, SeisPlot, SeisProcessing
using HCDSP, PyPlot, LinearAlgebra

# same as collect(ot:dt:tmax) or collect(range(start::T,stop::T,step::T=))
ot = 0.0; dt = 0.004; nt = 301; tmax = ot + (nt-1)*dt
taxis = range(ot,tmax,length=nt) |> collect; 

# set spatial dimension
ox = -1000.0; dx = 20; nx = 101; xmax = ox + (nx-1)*dx;
xaxis = range(ox,xmax,length=nx) |> collect;

# zero-offset travel-times
tau = [0.2, 0.6, 0.9]; 

# rms velocities
vel = [1500, 2000, 3000];

# apex-shifts (for dipping layers)
apex = zero(vel);

# events amplitudes
amp = [10, -1, 0.1];

# wavelet
wav_flag = "ricker";

# central freq
f0 = [20.0];

d = SeisHypEvents(f0 = f0,
                wavelet =  wav_flag,
                amp = amp,
                apex = apex,
                vel = vel,
                tau =  tau,
                ot = ot,
                ox = ox,
                nx = nx,
                nt = nt );
SeisPlotTX(d,pclip=95); gcf()