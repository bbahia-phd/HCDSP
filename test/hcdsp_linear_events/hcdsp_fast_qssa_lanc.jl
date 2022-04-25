cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using PyPlot
using SeisMain, SeisPlot
using HCDSP

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=10, ox2=0.0, dx2=20.0, nx2=10, ox3=0.0, dx3=10.0,
            nx3=10, ox4=0.0, dx4=10.0, nx4=10, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);
dnx = SeisAddNoise(dzx, -1.5, db=true, L=5);

Q = quaternion(dnx,dnx,dnx);

nin = size(Q); npad = nextpow.(2,nin);

#INF = PadOp(Q,nin=nin,npad=npad,flag="fwd");

OUTF = zero(Q);

INF = qfft(Q,μ0,"left",1);

fmin = 0.0; fmax = 100.0; dt = 0.004; k = 10;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( nin[2:end] ) ;

iω = 30;

inf = INF[iω,indx];

outf = Op(inf,30);

# dimensions of array
dims = size(inf)    

# compute L, K
L = floor.(Int64,dims ./ 2) .+ 1;
K = dims .- L .+ 1;

H = HCDSP.build_hankel_matrix(inf);

fwd(x) = qmbh_multiply(inf, x, side="left", flag="fwd");
v = zeros(eltype(inf),prod(K))
for i in 1:prod(K)
    v[i] = randq()
end
norm(fwd(v) .- H*v,2).^2

adj(x) = qmbh_multiply(inf, x, side="left", flag="adj");
t = zeros(eltype(inf),prod(L))
for i in 1:prod(L)
    t[i] = randq()
end
norm(adj(t) .- H'*t,2).^2

#####
function Op(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

    fwd(x) = mbh_multiply(IN,x,flag="fwd");
    adj(x) = mbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, Bk, V = lanbpro(A,k,m=L1*L2*L3*L4,n=K1*K2*K3*K4,qflag=true)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += anti_diagonal_summation(Ub[:,i],V[:,i],(L1,L2,L3,L4),(K1,K2,K3,K4));
    end

    count = HCDSP.count_copy_times([10,10,10,10])

    return out ./ count
end

