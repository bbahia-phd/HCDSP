cd(joinpath(homedir(),"projects"))
pwd()

using Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP/"))
Pkg.status()

using Revise

using FFTW
using PyPlot
using SeisMain, SeisPlot
using HCDSP, IterativeMethods

params_zx = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[1.0,-1.0], f0=20.0)
dzx = SeisLinearEvents(; params_zx...);
dnx = SeisAddNoise(dzx, 0.5, db=true, L=9);

params_zy = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[-1.0,1.0], f0=20.0)
dzy = SeisLinearEvents(; params_zy...);
dny = SeisAddNoise(dzy, 1.0, db=true, L=9);

params_zz = (ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
            nx1=50, ox2=0.0, dx2=10.0, nx2=50, ox3=0.0, dx3=10.0,
            nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
            p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
            amp=[-1.0,-1.0], f0=20.0)
dzz = SeisLinearEvents(; params_zz...);
dnz = SeisAddNoise(dzz, 3.0, db=true, L=9);

function core_sing_vals_dist(;ms=1.0)
    # Define time domain quaternion
    Q = quaternion(dzx,dzy,dzz);

    nin = size(Q); npad = nextpow.(2,nin);

    INX = PadOp(dzx,nin=nin,npad=npad,flag="fwd");
    INX = fft(INX,1);

    INY = PadOp(dzy,nin=nin,npad=npad,flag="fwd");
    INY = fft(INY,1);

    INZ = PadOp(dzz,nin=nin,npad=npad,flag="fwd");
    INZ = fft(INZ,1);

    INF = PadOp(Q,nin=nin,npad=npad,flag="fwd");
    INF = qfft(Q,μ0,"left",1);

    fmin = 0.0; fmax = 100.0; dt = 0.004; k = 4;

    # Freq range
    ω_range = freq_indexes(fmin, fmax, dt, npad[1])

    # Spatial indexes
    indx = CartesianIndices( nin[2:end] ) ;

    iω = 50;

    inx = INX[iω,indx];
    iny = INY[iω,indx];
    inz = INZ[iω,indx];

    inf = INF[iω,indx];

    # dimensions of array
    dims = size(inf)    

    # compute L, K
    L = floor.(Int64,dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    # trajectory matrix
    Hx = HCDSP.build_hankel_matrix(inx);
    Hy = HCDSP.build_hankel_matrix(iny);
    Hz = HCDSP.build_hankel_matrix(inz);

    Hq = HCDSP.build_hankel_matrix(inf);

    # Augmented trajectory
    Ha = vcat(Hq,invi.(Hq),invj.(Hq),invk.(Hq));

    # Complex adjoints
    Hc  = fwd_complex_adjoint(Hq);
    Hca = fwd_complex_adjoint(Ha);

    # Complex SVD
    Ux,Σx,Vx = LinearAlgebra.svd(Hx);
    Uy,Σy,Vy = LinearAlgebra.svd(Hy);
    Uz,Σz,Vz = LinearAlgebra.svd(Hz);

    Uq,Σq,Vq = LinearAlgebra.svd(Hc);
    Ua,Σa,Va = LinearAlgebra.svd(Hca);

    σq_plot = Σq[1:2:end];
    σa_plot = Σa[1:2:end];
    sing_axis = collect(1:length(σq_plot));

    fignum = "hcdsp_sing_vals"

    close("all"); clf()
    figure(fignum,figsize=(5,3))

    subplot(1,2,1)
    plot(sing_axis,σq_plot,marker="o",linestyle="",color="k",markersize=ms,label="Q");
    plot(sing_axis,σa_plot,marker="o",linestyle="",color="b",markersize=ms,label="A");

    plot(sing_axis,Σx,marker="3",linestyle="",color="k",markersize=ms,label="X");
    plot(sing_axis,Σy,marker="4",linestyle="",color="k",markersize=ms,label="Y");
    plot(sing_axis,Σz,marker="1",linestyle="",color="k",markersize=ms,label="Z");

    title("(a)",fontsize=10)
    xlabel("Index")
    ylabel("Singular values")
    legend(fontsize=10)

    # Define time domain (noisy) quaternion
    Q = quaternion(dnx,dny,dnz);

    INX = PadOp(dzx,nin=nin,npad=npad,flag="fwd");
    INX = fft(INX,1);

    INY = PadOp(dzy,nin=nin,npad=npad,flag="fwd");
    INY = fft(INY,1);

    INZ = PadOp(dzz,nin=nin,npad=npad,flag="fwd");
    INZ = fft(INZ,1);

    INF = PadOp(Q,nin=nin,npad=npad,flag="fwd");
    INF = qfft(Q,μ0,"left",1);

    fmin = 0.0; fmax = 100.0; dt = 0.004; k = 4;

    # Freq range
    ω_range = freq_indexes(fmin, fmax, dt, npad[1])

    # Spatial indexes
    indx = CartesianIndices( nin[2:end] ) ;

    iω = 50;

    inx = INX[iω,indx];
    iny = INY[iω,indx];
    inz = INZ[iω,indx];

    inf = INF[iω,indx];

    # trajectory matrix
    Hx = HCDSP.build_hankel_matrix(inx);
    Hy = HCDSP.build_hankel_matrix(iny);
    Hz = HCDSP.build_hankel_matrix(inz);

    Hq = HCDSP.build_hankel_matrix(inf);

    # Augmented trajectory
    Ha = vcat(Hq,invi.(Hq),invj.(Hq),invk.(Hq));

    # Complex adjoints
    Hc  = fwd_complex_adjoint(Hq);
    Hca = fwd_complex_adjoint(Ha);

    # Complex SVD
    Ux,Σx,Vx = LinearAlgebra.svd(Hx);
    Uy,Σy,Vy = LinearAlgebra.svd(Hy);
    Uz,Σz,Vz = LinearAlgebra.svd(Hz);

    Uq,Σq,Vq = LinearAlgebra.svd(Hc);
    Ua,Σa,Va = LinearAlgebra.svd(Hca);

    σq_plot = Σq[1:2:end];
    σa_plot = Σa[1:2:end];
    sing_axis = collect(1:length(σq_plot));

    subplot(1,2,2)
    plot(sing_axis,σq_plot,marker="o",linestyle="",color="k",markersize=ms,label="Q");
    plot(sing_axis,σa_plot,marker="o",linestyle="",color="b",markersize=ms,label="A");

    plot(sing_axis,Σx,marker="3",linestyle="",color="k",markersize=ms,label="X");
    plot(sing_axis,Σy,marker="4",linestyle="",color="k",markersize=ms,label="Y");
    plot(sing_axis,Σz,marker="1",linestyle="",color="k",markersize=ms,label="Z");

    title("(b)",fontsize=10)
    xlabel("Index")
    ylabel("Singular values")
    legend(fontsize=10)

    tight_layout()
 
    gcf()
end