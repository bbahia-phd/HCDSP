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

params_zx = (ot=0.0, dt=0.004, nt=100, ox1=0.0, dx1=10.0,
            nx1=20, ox2=0.0, dx2=20.0, nx2=20, ox3=0.0, dx3=10.0,
            nx3=20, ox4=0.0, dx4=10.0, nx4=20, tau=[0.1,0.25],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);
dnx = SeisAddNoise(dzx, -1.5, db=true, L=5);

nin = size(dnx); npad = nextpow.(2,nin);

INF = complex.(PadOp(dnx,nin=nin,npad=npad,flag="fwd"));
OUTF = zero(INF);

fft!(INF,1);

fmin = 0.0; fmax = 100.0; dt = 0.004; k = 10;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

for iω in ω_range
    OUTF[iω,indx] .= Op(INF[iω,indx],k...)
end

# Symmetries and ifft
conj_symmetry!(OUTF)

# Truncate
OUT = PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );

j = 10;
close("all"); clf();
SeisPlotTX([dzx[:,:,j,j,j] dnx[:,:,j,j,j] OUT[:,:,j,j,j]],fignum="panel",style="wiggles",wbox=10,hbox=5);gcf()

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

    U, Bk, V = lanbpro(A,k,m=L1*L2*L3*L4,n=K1*K2*K3*K4)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += anti_diagonal_summation2(Ub[:,i],V[:,i],(L1,L2,L3,L4),(K1,K2,K3,K4));
    end

    count = HCDSP.count_copy_times([32,32,32,32])

    return out ./ count
end