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

#####
function Op(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = mbh_multiply(IN,x,flag="fwd");
    adj(x) = mbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, Bk, V = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K))

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += anti_diagonal_summation(Ub[:,i],V[:,i],L,K);
    end

    count = HCDSP.count_copy_times([10,10,10,10])

    return out ./ count
end

drx = fx_process(dnx, 0.004, 0, 100, Op, (10,false)...);


#####
function QOp(IN,k,qflag=true)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = qmbh_multiply(IN,x,flag="fwd");
    adj(x) = qmbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, Bk, V = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K),qflag=qflag)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += anti_diagonal_summation(Ub[:,i],V[:,i],L,K);
    end

    count = HCDSP.count_copy_times([10,10,10,10])

    return out ./ count
end

# Quaternion denoising
Q = quaternion(dnx,dnx,dnx);
QQ = fx_process(Q, 0.004, 0, 100, QOp,(10,true)...);
dqx = imagi.(QQ);

# prediction_quality
Rx = quality(drx,dzx)

Qx = quality(dqx,dzx)

# get slice j
j = 5;

close("all"); clf();
SeisPlotTX([dzx[:,:,j,j,j] dnx[:,:,j,j,j] drx[:,:,j,j,j] dqx[:,:,j,j,j]],fignum="panel",style="wiggles",wbox=10);
gcf()