cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using HCDSP
using PyPlot
using SeisMain, SeisPlot

function LANCSSAOp(IN,k)

    # Hankelize
    H = HankelOp(IN);
    
    # Rank reduction
    U, Bk, V = HCDSP.lanbpro(H,k)
    
    # Averaging
    OUT = AveragingOp(U*U'*H,size(IN))

    return OUT
end

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);

nin = size(dzx); npad = nextpow.(2,nin);

INC = complex.(PadOp(dzx,nin=nin,npad=npad,flag="fwd"));

fft!(INC,1);

fmin = 0.0; fmax = 100.0; dt = 0.004;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

iω = 15;

dc = copy(INC[iω,indx]);

Op1(d,k) = SVDSSAOp(d,k);
Op2(d,k) = rQROp(d,k);
Op3(d,k) = LANCSSAOp(d,k);

K = 2:2:50;
kmax = length(K)
o1  = zeros(eltype(dc),size(dc)...,kmax);
o2  = zeros(eltype(dc),size(dc)...,kmax);
o3  = zeros(eltype(dc),size(dc)...,kmax);
tmp = zeros(eltype(dc),size(dc)...);

rmax = 20;

r1 = zeros(rmax,kmax);
r2 = zeros(rmax,kmax);
r3 = zeros(rmax,kmax);

for k in 1:kmax
    kk = K[k]

    for r in 1:rmax
        dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
        INF = complex.(PadOp(dnx,nin=nin,npad=npad,flag="fwd"));
        fft!(INF,1);

        d  = copy(INF[iω,indx]);

        tmp .= Op1(d,kk);
        r1[r,k] = quality(tmp,dc);
        #o1[:,:,kk] .= tmp;

        tmp .= Op2(d,kk);
        r2[r,k] = quality(tmp,dc);
        #o2[:,:,kk] .= tmp;

        tmp .= Op3(d,kk);
        r3[r,k] = quality(tmp,dc);
        #o3[:,:,kk] .= tmp;
    end
    @show [kk  sum(r1[:,k],dims=1) sum(r2[:,k],dims=1) sum(r3[:,k],dims=1)]
end

# Average
r1r = sum(r1,dims=1)' ./ rmax;
r2r = sum(r2,dims=1)' ./ rmax;
r3r = sum(r3,dims=1)' ./ rmax;