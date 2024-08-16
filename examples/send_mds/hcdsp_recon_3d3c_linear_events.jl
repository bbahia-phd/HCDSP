
# Change this path to your HCDSP path
cd("/lustre03/vol0/4638ns/projects/HCDSP/")
pwd()

using Distributed
addprocs(20)

@everywhere using Pkg

# Change this path to your HCDSP path
# You will need to  load dependencies by using Pkg.instantiate().
# This step might need some troubleshooting, please let me know.
@everywhere Pkg.activate("/lustre03/vol0/4638ns/projects/HCDSP/")

@everywhere using Revise            # To avoid reinitializing julia
@everywhere using LinearAlgebra     
@everywhere using FFTW
@everywhere using Random

@everywhere using HCDSP

# Note that I am not using SeisPlot.jl because there is a dependency conflict.
# Not sure what has changed, but you might just use SeisMakie.jl.
# Still, I copied all SeisPlot.jl functions into HCDSP so that you can just use SeisPlotTX out the box here.
using PyPlot
using SeisMain
using HDF5

# Creating pure mode synthetics
function get_mode_data(;nx1=40,nx2=40,nx3=1,nx4=1)

    params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.1],
    p1=[0.0001],p2=[-0.0001],p3=[0.0002],p4=[-0.0002],
    amp=[1.0], f0=20.0)
    p = SeisLinearEvents(; params_zx...);

    params_zy = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.25],
    p1=[-0.0003],p2=[0.0003],p3=[-0.0001],p4=[0.0002],
    amp=[-1.0], f0=20.0)
    sv = SeisLinearEvents(; params_zy...);

    params_zz = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.3],
    p1=[-0.0002],p2=[0.0001],p3=[-0.0003],p4=[0.0001],
    amp=[-1.0], f0=20.0)
    sh = SeisLinearEvents(; params_zz...);

    return (p,sv,sh)

end

# Unmixing components (assumes known "angles of incidence". Perhaps a better example would be good too.)
function unmix(p,sv,sh)

    A = inv([0.75 0.15 0.1; 0.15 0.75 0.1; 0.1 0.15 0.75]);
    
    o1,o2,o3 = similar(p),similar(p),similar(p)
    for i in eachindex(p)
        tmp = A*[p[i]; sv[i]; sh[i]]
        o1[i] = tmp[1];
        o2[i] = tmp[2];
        o3[i] = tmp[3];
    end

    return o1,o2,o3
end

# Mixing components (assumes known "angles of incidence". Perhaps a better example would be good too.)
function mix(p,sv,sh)

    o1,o2,o3 = similar(p),similar(p),similar(p)
    for i in eachindex(p)
        o1[i] = 0.75*p[i] + 0.15*sv[i] + 0.1*sh[i];
        o2[i] = 0.15*p[i] + 0.75*sv[i] + 0.1*sh[i];
        o3[i] = 0.1*p[i]  + 0.15*sv[i] + 0.75*sh[i];
    end

    return o1,o2,o3
end

# clean & pure seismic modes
p,sv,sh = get_mode_data(nx1=20,nx2=20,nx3=20,nx4=20);

# mixed observed displacements
# might be considered a mid-to-far offset assumption (nears would have close-to-vertical rays)
dzz,dzy,dzx = mix(p,sv,sh);

# fx-process setup
fmin = 0.0; fmax = 60.0; dt = 0.004;
@everywhere α = 0.5;

# Define operator to act on a frequency slice d
@everywhere imp_ssa(d,k)   = HCDSP.imputation_op(d,HCDSP.fast_ssa_lanc,  (k); iter=100, α = α)
@everywhere imp_qssa(d,k)  = HCDSP.imputation_op(d,HCDSP.fast_qssa_lanc, (k); iter=100, α = α)
@everywhere imp_aqssa(d,k) = HCDSP.imputation_op(d,HCDSP.fast_aqssa_lanc,(k); iter=100, α = α)

# SNRs
snrx,snry,snrz=0.8,1.0,1.2;

# Add noise
dnx = SeisAddNoise(dzx, snrx, db=true, L=3);
dny = SeisAddNoise(dzy, snry, db=true, L=3);
dnz = SeisAddNoise(dzz, snrz, db=true, L=3);

# Temporary Quaternion
Qt = quaternion(dnx,dny,dnz);

# decimations to test
perc = 90;

# Missing traces
Qt .= decimate_traces(Qt,perc);

# ranks to test
k  = 10;
ka = 12;

# Call pmap fx_process with Q imputation (runs on both sides of the spectra)
Qo = pmap_fx_process(Qt,dt,fmin,fmax,imp_qssa,(k));
qx = quality(imagi.(Qo),dzx)
qy = quality(imagj.(Qo),dzy)
qz = quality(imagk.(Qo),dzz)

# Call pmap fx_process with Q imputation (runs on both sides of the spectra)
Qa = pmap_fx_process(Qt,dt,fmin,fmax,imp_aqssa,(ka));
aqx = quality(imagi.(Qa),dzx)
aqy = quality(imagj.(Qa),dzy)
aqz = quality(imagk.(Qa),dzz)

# Component-wise processing (runs on single side of the spectra)
Xo = pmap_fx_process(imagi.(Qt),dt,fmin,fmax,imp_ssa,(k));
Yo = pmap_fx_process(imagj.(Qt),dt,fmin,fmax,imp_ssa,(k));
Zo = pmap_fx_process(imagk.(Qt),dt,fmin,fmax,imp_ssa,(k));

# Get quality
rx = quality(Xo,dzx)
ry = quality(Yo,dzy)
rz = quality(Zo,dzz)

n=15;
clf();close("all")
SeisPlotTX(
    [dzx[:,:,n] imagi.(Qt)[:,:,n] Xo[:,:,n] imagi.(Qo)[:,:,n] imagi.(Qa)[:,:,n] (dzx .- Xo)[:,:,n] (dzx .- imagi.(Qo))[:,:,n] (dzx .- imagi.(Qa))[:,:,n]], wbox=20,  hbox=4, cmap="gray",xcur=2.0);
gcf()
# PyPlot.savefig(joinpath(homedir(),"julia/compare_x.jpg"));

clf();close("all")
SeisPlotTX(
    [dzy[:,:,n] imagj.(Qt)[:,:,n] Yo[:,:,n] imagj.(Qo)[:,:,n] imagj.(Qa)[:,:,n] (dzy .- Yo)[:,:,n] (dzy .- imagj.(Qo))[:,:,n] (dzy .- imagj.(Qa))[:,:,n]], wbox=20,  hbox=4, cmap="gray",xcur=2.0);
gcf()
# PyPlot.savefig(joinpath(homedir(),"julia/compare_y.jpg"));

clf();close("all")
SeisPlotTX(
    [dzz[:,:,n] imagk.(Qt)[:,:,n] Zo[:,:,n] imagk.(Qo)[:,:,n] imagk.(Qa)[:,:,n] (dzz .- Zo)[:,:,n] (dzz .- imagk.(Qo))[:,:,n] (dzz .- imagk.(Qa))[:,:,n]], wbox=20,  hbox=4, cmap="gray",xcur=2.0);
gcf()
# PyPlot.savefig(joinpath(homedir(),"julia/compare_z.jpg"));