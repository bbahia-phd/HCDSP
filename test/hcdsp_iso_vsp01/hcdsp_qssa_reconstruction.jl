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
    H = vcat( H,
            invi.(H),
            invj.(H),
            invk.(H) );

    # Rank reduction
    U, Bk, V = HCDSP.lanbpro(H,k)
    
    # Averaging
    OUT = AveragingOp(U*U'*H,size(IN))

    return OUT
end

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

    A = inv([0.75 0.15 0.05; 0.15 0.75 0.05; 0.05 0.15 0.75]);
    
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
        o1[i] = 0.75*p[i] + 0.15*sv[i] + 0.05*sh[i];
        o2[i] = 0.15*p[i] + 0.75*sv[i] + 0.05*sh[i];
        o3[i] = 0.05*p[i] + 0.15*sv[i] + 0.75*sh[i];
    end

    return o1,o2,o3
end

# clean & pure seismic modes
p,sv,sh = get_mode_data();

# mixed observed displacements
dzz,dzy,dzx = mix(p,sv,sh);

# noisy displacements
dnx = SeisAddNoise(dzx, -1.0, db=true, L=3);
dny = SeisAddNoise(dzy, -1.0, db=true, L=3);
dnz = SeisAddNoise(dzz, -1.0, db=true, L=3);

# FSSA
nin = size(dnx); npad = nextpow.(2,nin); dt = 0.004; fmin = 0.0; fmax = 50; k = 6;

# Quaternion denoising
Qc = Quaternion.(dzx,dzy,dzz);
Qn = Quaternion.(dnx,dny,dnz);

QC = PadOp(Qc,nin=nin,npad=npad,flag="fwd");
QN = PadOp(Qn,nin=nin,npad=npad,flag="fwd");

# qfts
QC = fft(QC,1);
QN = fft(QN,1);

# 
fmin = 0.0; fmax = 100.0; dt = 0.004; k = 10;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

iω = 10;

d  = copy(QN[iω,indx]);
dc = copy(QC[iω,indx]);

Op1(d,k) = SVDSSAOp(d,k);
Op2(d,k) = rQROp(d,k);
Op3(d,k) = LANCSSAOp(d,k)

kmax = 25;
o1  = zeros(eltype(d),size(d)...,kmax);
o2  = zeros(eltype(d),size(d)...,kmax);
o3  = zeros(eltype(d),size(d)...,kmax);
tmp = zeros(eltype(d),size(d)...);

rmax = 1;

r1 = zeros(rmax,kmax);
r2 = zeros(rmax,kmax);
r3 = zeros(rmax,kmax);

for kk in 1:kmax

    for r in 1:rmax

        # noisy displacements
        dnx = SeisAddNoise(dzx, -1.0, db=true, L=3);
        dny = SeisAddNoise(dzy, -1.0, db=true, L=3);
        dnz = SeisAddNoise(dzz, -1.0, db=true, L=3);

        # Quaternion denoising
        Qn = Quaternion.(dnx,dny,dnz);
        QN = PadOp(Qn,nin=nin,npad=npad,flag="fwd");
        QN .= fft(QN,1);

        iω = 10;

        d  = copy(QN[iω,indx]);

        tmp .= Op1(d,kk);
        r1[r,kk] = quality(tmp,dc);
        #o1[:,:,kk] .= tmp;

        tmp .= Op2(d,kk);
        r2[r,kk] = quality(tmp,dc);
        #o2[:,:,kk] .= tmp;

        tmp .= Op3(d,kk);
        r3[r,kk] = quality(tmp,dc);
        #o3[:,:,kk] .= tmp;
    end

    @show [kk sum(r1[:,kk],dims=1) sum(r2[:,kk],dims=1) sum(r3[:,kk],dims=1)]
end

r1q = sum(r1,dims=1)' ./ rmax;
r2q = sum(r2,dims=1)' ./ rmax;
r3q = sum(r3,dims=1)' ./ rmax;