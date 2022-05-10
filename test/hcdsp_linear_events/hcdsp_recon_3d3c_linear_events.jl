pwd()

dev_dir="/dev/Breno_GOM/projects/";

cd(dev_dir)
#cd(joinpath(homedir(),"projects"))

pwd()

using Pkg
Pkg.activate(joinpath(dev_dir,"HCDSP/"))
Pkg.status()

using Revise

using FFTW
using PyPlot
using LinearAlgebra
using StatsBase,Statistics
using SeisMain, SeisPlot
using HCDSP


function get_mode_data(;nx1=40,nx2=40,nx3=1,nx4=1)

    params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.1],
    p1=[0.0001],p2=[0.0],p3=[0.0],p4=[0.0],
    amp=[1.0], f0=20.0)
    p = SeisLinearEvents(; params_zx...);

    params_zy = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.25],
    p1=[-0.0003],p2=[0.0],p3=[0.0],p4=[0.0],
    amp=[-1.0], f0=20.0)
    sv = SeisLinearEvents(; params_zy...);

    params_zz = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
    nx1=nx1, ox2=0.0, dx2=10.0, nx2=nx2, ox3=0.0, dx3=10.0,
    nx3=nx3, ox4=0.0, dx4=10.0, nx4=nx4, tau=[0.3],
    p1=[-0.0002],p2=[0.0],p3=[0.0],p4=[0.0],
    amp=[-1.0], f0=20.0)
    sh = SeisLinearEvents(; params_zz...);

    return (p,sv,sh)

end

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
p,sv,sh = get_mode_data();

# mixed observed displacements
dzz,dzy,dzx = mix(p,sv,sh);

fmin = 0.0; fmax = 80.0; dt = 0.004;

# Define operator to act on a frequency slice d
imp_ssa(d,k)   = HCDSP.imputation_op(d,HCDSP.fast_ssa_lanc,  (k); iter=10)
imp_qssa(d,k)  = HCDSP.imputation_op(d,HCDSP.fast_qssa_lanc, (k); iter=10)
imp_aqssa(d,k) = HCDSP.imputation_op(d,HCDSP.fast_aqssa_lanc,(k); iter=10)

# SNRs
snrx,snry,snrz=0.8,1.0,1.2;

# Add noise
dnx = SeisAddNoise(dzx, snrx, db=true, L=3);
dny = SeisAddNoise(dzy, snry, db=true, L=3);
dnz = SeisAddNoise(dzz, snrz, db=true, L=3);

# Temporary Quaternion
Qt = quaternion(dnx,dny,dnz);

# decimations to test
perc = 50;

# Missing traces
Qt .= decimate_traces(Qt,perc);

# ranks to test
k = 8;

# Component-wise processing
Xo = fx_process(imagi.(Qt),dt,fmin,fmax,imp_ssa,(k));
Yo = fx_process(imagj.(Qt),dt,fmin,fmax,imp_ssa,(k));
Zo = fx_process(imagk.(Qt),dt,fmin,fmax,imp_ssa,(k));

# Call fx_process with Q imputation
Qo = fx_process(Qt,dt,fmin,fmax,imp_qssa,(k));
Qa = fx_process(Qt,dt,fmin,fmax,imp_aqssa,(12));

# Get quality
rx = quality(Xo,dzx)
ry = quality(Yo,dzy)
rz = quality(Zo,dzz)

qx = quality(imagi.(Qo),dzx)
qy = quality(imagj.(Qo),dzy)
qz = quality(imagk.(Qo),dzz)

aqx = quality(imagi.(Qa),dzx)
aqy = quality(imagj.(Qa),dzy)
aqz = quality(imagk.(Qa),dzz)

n=10;
clf();close("all")
SeisPlotTX(
    [dzx[:,:,n] imagi.(Qt)[:,:,n] Xo[:,:,n] imagi.(Qo)[:,:,n] imagi.(Qa)[:,:,n] (dzx .- Xo)[:,:,n] (dzx .- imagi.(Qo))[:,:,n] (dzx .- imagi.(Qa))[:,:,n]], wbox=8,  hbox=4, cmap="gray",xcur=2.0);
gcf()

clf();close("all")
SeisPlotTX(
    [dzy[:,:,n] imagj.(Qt)[:,:,n] Yo[:,:,n] imagj.(Qo)[:,:,n] imagj.(Qa)[:,:,n] (dzy .- Yo)[:,:,n] (dzy .- imagj.(Qo))[:,:,n] (dzy .- imagj.(Qa))[:,:,n]], wbox=8,  hbox=4, cmap="gray",xcur=2.0);
gcf()

clf();close("all")
SeisPlotTX(
    [dzz[:,:,n] imagk.(Qt)[:,:,n] Zo[:,:,n] imagk.(Qo)[:,:,n] imagk.(Qa)[:,:,n] (dzz .- Zo)[:,:,n] (dzz .- imagk.(Qo))[:,:,n] (dzz .- imagk.(Qa))[:,:,n]], wbox=8,  hbox=4, cmap="gray",xcur=2.0);
gcf()


file_path = "./HCDSP/data/hcdsp_recon_3d3c_linear_events"

## Write clean data to file
file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_zz.bin")
read_write(file,"w";n=size(dzz),input=dzz,T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_zy.bin")
read_write(file,"w";n=size(dzy),input=dzy,T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_zx.bin")
read_write(file,"w";n=size(dzx),input=dzx,T=Float32)

## Write Input data to file
file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_noisy_zz.bin")
read_write(file,"w";n=size(dzz),input=imagk.(Qt),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_noisy_zy.bin")
read_write(file,"w";n=size(dzy),input=imagj.(Qt),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_noisy_zx.bin")
read_write(file,"w";n=size(dzx),input=imagi.(Qt),T=Float32)

## Write SSA output to file
file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_ssa_zz.bin")
read_write(file,"w";n=size(dzz),input=Zo,T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_ssa_zy.bin")
read_write(file,"w";n=size(dzy),input=Yo,T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_ssa_zx.bin")
read_write(file,"w";n=size(dzx),input=Xo,T=Float32)

## Write QSSA output to file
file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_qssa_zz.bin")
read_write(file,"w";n=size(dzz),input=imagk.(Qo),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_qssa_zy.bin")
read_write(file,"w";n=size(dzy),input=imagj.(Qo),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_qssa_zx.bin")
read_write(file,"w";n=size(dzx),input=imagi.(Qo),T=Float32)

## Write AQSSA output to file
file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_aqssa_zz.bin")
read_write(file,"w";n=size(dzz),input=imagk.(Qa),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_aqssa_zy.bin")
read_write(file,"w";n=size(dzy),input=imagj.(Qa),T=Float32)

file = joinpath(file_path,"hcdsp_recon_3d3c_linear_events_aqssa_zx.bin")
read_write(file,"w";n=size(dzx),input=imagi.(Qa),T=Float32)