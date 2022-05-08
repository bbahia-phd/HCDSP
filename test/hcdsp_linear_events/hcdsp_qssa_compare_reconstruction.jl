cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using HCDSP
using PyPlot
using LinearAlgebra
using StatsBase,Statistics
using SeisMain, SeisPlot

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[0.9,-0.3], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);

params_zy = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[-0.1,0.7], f0=20.0)
dzy = SeisLinearEvents(; params_zy...);

params_zz = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[0.1,-0.7], f0=20.0)
dzz = SeisLinearEvents(; params_zz...);

Q = quaternion(dzx,dzy,dzz);

nin = size(dzx); npad = nin;

INC = PadOp(Q,nin=nin,npad=npad,flag="fwd");
INC .= fft(INC,1);

fmin = 0.0; fmax = 50.0; dt = 0.004;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

iω = 15;

dc = copy(INC[iω,indx]);

# These act on a frequency slice d
imp_ssa(d,T,k) = HCDSP.imputation_op(d,T,SVDSSAOp,(k); iter=10)
imp_rqr(d,T,k) = HCDSP.imputation_op(d,T,rQROp,(k);    iter=10)
imp_lan(d,T,k) = HCDSP.imputation_op(d,T,LANCSSAOp,(k);iter=10)

# ranks to test
K = 1:2:50;
kmax = length(K);

# decimations to test
percs = 10:10:90;
pmax = length(percs);

# for outputs
tmp = zeros(eltype(dc),size(dc)...);

# number of realizations
rmax = 20;

r1x = zeros(rmax,pmax,kmax);
r1y = zeros(rmax,pmax,kmax);
r1z = zeros(rmax,pmax,kmax);

r2x = zeros(rmax,pmax,kmax);
r2y = zeros(rmax,pmax,kmax);
r2z = zeros(rmax,pmax,kmax);

r3x = zeros(rmax,pmax,kmax);
r3y = zeros(rmax,pmax,kmax);
r3z = zeros(rmax,pmax,kmax);

# Add noise
snrx,snry,snrz=0.8,1.0,1.2;

for k in 1:kmax
    kk = K[k]

    for p in 1:pmax
        perc = percs[p]

        for r in 1:rmax
            # Noise
            dnx = SeisAddNoise(dzx, snrx, db=true, L=3);
            dny = SeisAddNoise(dzy, snry, db=true, L=3);
            dnz = SeisAddNoise(dzz, snrz, db=true, L=3);            

            # Temporary Quaternion
            Qt = quaternion(dnx,dny,dnz);
            
            # Missing traces
            Qt .= decimate_traces(Qt,perc);
            
            # FFT
            INF = PadOp(Qt,nin=nin,npad=npad,flag="fwd");
            INF .= fft(INF,1);

            # Extract slice
            d  = copy(INF[iω,indx]);

            # Get sampling
            T = SamplingOp(d);
    
            tmp .= imp_ssa(d,T,kk);
            r1x[r,p,k] = quality(imagi.(tmp),imagi.(dc));
            r1y[r,p,k] = quality(imagj.(tmp),imagj.(dc));
            r1z[r,p,k] = quality(imagk.(tmp),imagk.(dc));
    
            tmp .= imp_rqr(d,T,kk);
            r2x[r,p,k] = quality(imagi.(tmp),imagi.(dc));
            r2y[r,p,k] = quality(imagj.(tmp),imagj.(dc));
            r2z[r,p,k] = quality(imagk.(tmp),imagk.(dc));
    
            tmp .= imp_lan(d,T,kk);
            r3x[r,p,k] = quality(imagi.(tmp),imagi.(dc));
            r3y[r,p,k] = quality(imagj.(tmp),imagj.(dc));
            r3z[r,p,k] = quality(imagk.(tmp),imagk.(dc));
        end       
        @show [kk perc mean(r1[:,p,k],dims=1) mean(r2[:,p,k],dims=1) mean(r3[:,p,k],dims=1)]    
    end
end

using HDF5

fname = joinpath(homedir(),"projects/HCDSP/data/hcdsp_qssa_compare_reconstruction.h5")
fid = h5open(fname, "w")

create_group(fid,"gains")
fid["gains"]["svd"] = r1;
fid["gains"]["rqr"] = r2;
fid["gains"]["lanc"] = r3;

close(fid)

#EOF