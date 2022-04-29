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
        o3[i] = 0.1*p[i] + 0.15*sv[i] + 0.75*sh[i];
    end

    return o1,o2,o3
end

# clean & pure seismic modes
p,sv,sh = get_mode_data();

# mixed observed displacements
dzz,dzy,dzx = mix(p,sv,sh);

# noisy displacements
dnx = SeisAddNoise(dzx, -2.5, db=true, L=3);
dny = SeisAddNoise(dzy, -1.0, db=true, L=3);
dnz = SeisAddNoise(dzz,  1.0, db=true, L=3);

# Quaternion denoising
Q = quaternion(dnx,dny,dnz);

dt = 0.004; fmin = 0.0; fmax = 100; k = 3; # best trial-and-error

QQ = fx_process(Q, dt, fmin, fmax, SVDSSAOp,(2k)...);
dqx = imagi.(QQ);dqy = imagj.(QQ); dqz = imagk.(QQ);
pq,svq,shq = unmix(dqz,dqy,dqx);

drx = fx_process(dnx, dt, fmin, fmax, SVDSSAOp, (k)...);
dry = fx_process(dny, dt, fmin, fmax, SVDSSAOp, (k)...);
drz = fx_process(dnz, dt, fmin, fmax, SVDSSAOp, (k)...);
pr,svr,shr = unmix(drz,dry,drx);

quality(dqz,dzz) 
quality(drz,dzz)

quality(dqy,dzy)
quality(dry,dzy)

quality(dqx,dzx)
quality(drx,dzx)

# get slice j
j = 5;

close("all"); clf();
SeisPlotTX([dzz[:,:,j] dnz[:,:,j] dqz[:,:,j] drz[:,:,j]],fignum="panel",style="overlay",wbox=12,hbox=5,xcur=5.0);
gcf()

close("all"); clf();
SeisPlotTX([dzy[:,:,j] dny[:,:,j] dqy[:,:,j] dry[:,:,j]],fignum="panel",style="overlay",wbox=12,hbox=5,xcur=5.0);
gcf()

close("all"); clf();
SeisPlotTX([dzx[:,:,j] dnx[:,:,j] dqx[:,:,j] drx[:,:,j]],fignum="panel",style="overlay",wbox=12,hbox=5,xcur=5.0);
gcf()

quality(pq,p) 
quality(pr,p)

quality(svq,sv)
quality(svr,sv)

quality(shq,sh)
quality(shr,sh)

close("all"); clf();
SeisPlotTX([p[:,:,j] pq[:,:,j] pr[:,:,j]],fignum="panel",style="overlay",wbox=10);
gcf()

close("all"); clf();
SeisPlotTX([sv[:,:,j] svq[:,:,j] svr[:,:,j]],fignum="panel",style="overlay",wbox=10,xcur=1.5,pclip=95);
gcf()

close("all"); clf();
SeisPlotTX([sh[:,:,j] shq[:,:,j] shr[:,:,j]],fignum="panel",style="overlay",wbox=10,xcur=1.5,pclip=95);
gcf()