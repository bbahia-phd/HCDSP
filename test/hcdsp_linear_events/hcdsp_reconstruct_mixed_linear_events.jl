# cd(joinpath(homedir(),"projects"))
# pwd()

# using Pkg
# Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
# Pkg.status()

# using Revise

# using FFTW
# using HCDSP
# using PyPlot
# using LinearAlgebra
# using StatsBase,Statistics
# using SeisMain, SeisPlot

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

# ranks to test
K = 1:2:50;
kmax = length(K);

# decimations to test
percs = 10:10:90;
pmax = length(percs);

# number of realizations
rmax = 10;

# reconstruction gains
rx = zeros(rmax,pmax,kmax);
ry = zeros(rmax,pmax,kmax);
rz = zeros(rmax,pmax,kmax);

qx = zeros(rmax,pmax,kmax);
qy = zeros(rmax,pmax,kmax);
qz = zeros(rmax,pmax,kmax);

aqx = zeros(rmax,pmax,kmax);
aqy = zeros(rmax,pmax,kmax);
aqz = zeros(rmax,pmax,kmax);

# SNRs
snrx,snry,snrz=0.8,1.0,1.2;

for k in 1:kmax
    # Get the rank
    kk = K[k]

    for p in 1:pmax

        # Get the decimation percentage
        perc = percs[p]

        # Start realizations
        for r in 1:rmax

            # Add noise
            dnx = SeisAddNoise(dzx, snrx, db=true, L=3);
            dny = SeisAddNoise(dzy, snry, db=true, L=3);
            dnz = SeisAddNoise(dzz, snrz, db=true, L=3);
            
            # Temporary Quaternion
            Qt = quaternion(dnx,dny,dnz);
            
            # Missing traces
            Qt .= decimate_traces(Qt,perc);

            # Component-wise processing
            Xo = fx_process(imagi.(Qt),dt,fmin,fmax,imp_ssa,(k))
            Yo = fx_process(imagj.(Qt),dt,fmin,fmax,imp_ssa,(k))
            Zo = fx_process(imagk.(Qt),dt,fmin,fmax,imp_ssa,(k))

            # Call fx_process with Q imputation
            Qo = fx_process(Qt,dt,fmin,fmax,imp_qssa,(k))
            Qa = fx_process(Qt,dt,fmin,fmax,imp_aqssa,(k))
            
            # Get quality
            rx[r,p,k] = quality(Xo,dzx);
            ry[r,p,k] = quality(Yo,dzy);
            rz[r,p,k] = quality(Zo,dzz);

            qx[r,p,k] = quality(imagi.(Qo),dzx);
            qy[r,p,k] = quality(imagj.(Qo),dzy);
            qz[r,p,k] = quality(imagk.(Qo),dzz);

            aqx[r,p,k] = quality(imagi.(Qa),dzx);
            aqy[r,p,k] = quality(imagj.(Qa),dzy);
            aqz[r,p,k] = quality(imagk.(Qa),dzz);
        end       
        @show [kk perc mean(rx[:,p,k],dims=1) mean(ry[:,p,k],dims=1) mean(rz[:,p,k],dims=1)]
        @show [kk perc mean(qx[:,p,k],dims=1) mean(qy[:,p,k],dims=1) mean(qz[:,p,k],dims=1)]                
        @show [kk perc mean(aqx[:,p,k],dims=1) mean(aqy[:,p,k],dims=1) mean(aqz[:,p,k],dims=1)]                
    end
end

using HDF5

file_path = "./HCDSP/data/hcdsp_recon_3d3c_linear_events"

fname = joinpath(file_path,"hcdsp_reconstruct_mixed_linear_events.h5")
fid = h5open(fname, "w")

create_group(fid,"gains/real")
fid["gains/real"]["x"] = rx;
fid["gains/real"]["y"] = ry;
fid["gains/real"]["z"] = rz;

create_group(fid,"gains/quater")
fid["gains/quater"]["x"] = qx;
fid["gains/quater"]["y"] = qy;
fid["gains/quater"]["z"] = qz;

create_group(fid,"gains/quater")
fid["gains/aquater"]["x"] = aqx;
fid["gains/aquater"]["y"] = aqy;
fid["gains/aquater"]["z"] = aqz;

close(fid)