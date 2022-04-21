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

# FSSA
nin = size(dnx); npad = nin; dt = 0.004; fmin = 0.0; fmax = 50; k = 6;

# Quaternion denoising
Qc = Quaternion.(dzx,dzy,dzz);
QC = PadOp(Qc,nin=nin,npad=npad,flag="fwd");

# qfts
QC = fft(QC,1);
QN = fft(QN,1);

# parameters 
fmin = 0.0; fmax = 100.0; dt = 0.004;

# Freq range
ω_range = freq_indexes(fmin, fmax, dt, npad[1])

# choose a frequency
iω = 10;

# Spatial indexes
indx = CartesianIndices( npad[2:end] ) ;

# get data for that frequency at every spatial index
dc = copy(QC[iω,indx]);

# Define the SSA Operators
Op1(d,k) = SVDSSAOp(d,k);
Op2(d,k) = rQROp(d,k);
Op3(d,k) = LANCSSAOp(d,k)

# Rank to test
K = 1:2:50;
kmax = length(K);

# Max nubmer of realizations
rmax = 10;

# a temp array
tmp = zeros(eltype(dc),size(dc)...);

# allocate
r1 = zeros(rmax,kmax);
r2 = zeros(rmax,kmax);
r3 = zeros(rmax,kmax);

for k in 1:kmax
    kk = K[k]    

    for r in 1:rmax
        # noisy displacements
        dnx = SeisAddNoise(dzx, -2.0, db=true, L=3);
        dny = SeisAddNoise(dzy, -2.0, db=true, L=3);
        dnz = SeisAddNoise(dzz, -2.0, db=true, L=3);

        # Quaternion denoising
        Qn = Quaternion.(dnx,dny,dnz);
        QN = PadOp(Qn,nin=nin,npad=npad,flag="fwd");
        QN .= fft(QN,1);

        # this was defined above
        iω = 10;

        # extract f slice
        d  = copy(QN[iω,indx]);

        # SVD
        tmp .= Op1(d,kk);
        r1[r,k] = quality(tmp,dc);

        # rQR
        tmp .= Op2(d,kk);
        r2[r,k] = quality(tmp,dc);

        # Lanczos
        tmp .= Op3(d,kk);
        r3[r,k] = quality(tmp,dc);
    end

    @show [kk sum(r1[:,k],dims=1) sum(r2[:,k],dims=1) sum(r3[:,k],dims=1)]
end

using HDF5

fname = "./hcdsp_qssa_compare_denoising"
fid = h5open(fname, "w")

create_group(fid,"gains")
fid["gains"]["svd"] = r1;
fid["gains"]["rqr"] = r2;
fid["gains"]["lanc"] = r3;

# close(fid)

#=
# Average
r1r = vec(mean(r1,dims=1));
r2r = vec(mean(r2,dims=1));
r3r = vec(mean(r3,dims=1));

std1 = vec(std(r1,dims=1));
std2 = vec(std(r2,dims=1));
std3 = vec(std(r3,dims=1));

close("all");
clf();

errorbar(K,r1r,yerr=std1,fmt="-o",label="SVD")
errorbar(K,r2r,yerr=std2,fmt="-o",label="rQR")
errorbar(K,r3r,yerr=std3,fmt="-o",label="Lanczos")
xlabel("'Rank' "*L"(k)",fontsize=15)
ylabel(L"R = 10\log{\left( \frac{ \parallel {\bf d}^{o} \parallel^2_F }{\parallel \hat{\bf d}_j - {\bf d}^{o} \parallel^2_F} \right)}",fontsize=15)
legend()

gcf()

close("all");
clf();
plot(K,r1r);
plot(K,r2r);
plot(K,r3r);
gcf()

close("all");clf()
fignum="time_domain_data";
figure(fignum,figsize=(9,5));

subplot(1,2,1);
SeisPlotTX(dzx[:,1:2:80,1],
        dy=dt,
        dx=10,
        cmap="gray",
        style="overlay",
        fignum=fignum,
        xcur=2,
        title="(a)",
        ylabel="Time (s)",
        xlabel="Offset (m)")

subplot(1,2,2);
SeisPlotTX(dnx[:,1:2:80,1],
        dy=dt,
        dx=10,
        cmap="gray",
        style="overlay",
        fignum=fignum,
        xcur=2,
        title="(b)",
        ylabel="Time (s)",
        xlabel="Offset (m)")
tight_layout();
gcf()
=#
