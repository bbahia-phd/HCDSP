pwd()

using Distributed
addprocs(19)

@everywhere dev_dir="/dev/Breno_GOM/projects/";

cd(dev_dir)
#cd(joinpath(homedir(),"projects"))

pwd()

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP/"))
@everywhere Pkg.status()

@everywhere using Revise

@everywhere using FFTW
@everywhere using PyPlot
@everywhere using LinearAlgebra
@everywhere using SeisMain, SeisPlot
@everywhere using HCDSP

@everywhere include(joinpath(dev_dir,"HCDSP/test/hcdsp_iso_vsp01/hcdsp_pmap_functions.jl"))

# data dir home
data_path = "/dev/Breno_GOM/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zx_crg1350.seis");
dzy,hzy,ext_zy = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zy_crg1350.seis");
dzz,hzz,ext_zz = SeisRead("./../iso_vsp01/sgy/iso_vsp01_zz_crg1350.seis");

dt      = Float64(ext_zz.d1);
nt      = Int64(ext_zz.n1);
ntr     = ext_zz.n2;
rz_init = 1350.0;  # m
rz_end  = 1850.01; # m
drz     = 16.667 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr      = length(rz_axis)
ns      = 205; # number of sources within lines
nsline  = 205; # number of source lines

@assert ns*nsline == ntr == size(dzz,2)

dzx = reshape(dzx,(nt,ns,nsline));
dzy = reshape(dzy,(nt,ns,nsline));
dzz = reshape(dzz,(nt,ns,nsline));


@everywhere begin
    
    # f-x params
    dt    = 0.012012;
    fmin  = 0;
    fmax  = 50;
    
    # patching params
    psize = (50,20,20);
    polap = (50,50,50);
    smin  = (1,1,1);
    smax  = (217,205,205);

    # rank-reduction params
    K = 1:2:10;
    iter = 20;
    α=0.5;

    # Define operator to act on a frequency slice d
    imp_ssa(d,k) = HCDSP.imputation_op(d,
                            HCDSP.fast_ssa_lanc,
                            (k); iter=iter,α=α)

    # Define operator to act on a frequency slice d
    imp_qssa(d,k) = HCDSP.imputation_op(d,
                            HCDSP.fast_qssa_lanc,
                            (k); iter=iter,α=α)

    # Define operator to act on a frequency slice d
    imp_aqssa(d,k) = HCDSP.imputation_op(d,
                            HCDSP.fast_aqssa_lanc,
                            (k); iter=iter,α=α)                            
  
end

# Add noise
snrx,snry,snrz=0.8,1.0,1.2;
dnx = SeisAddNoise(dzx, snrx, db=true,L=3);
dny = SeisAddNoise(dzy, snry, db=true,L=3);
dnz = SeisAddNoise(dzz, snrz, db=true,L=3);

# Decimate quaternion: all have same traces missing
perc = 50;
Qc = quaternion(dzx,dzy,dzz);
Qt = quaternion(dnx,dny,dnz);
Qt = decimate_traces(Qt,perc);

# split inputs
dx = imagi.(Qt);
dy = imagj.(Qt);
dz = imagk.(Qt);

Rx,Ry,Rz = zeros(length(K)),zeros(length(K)),zeros(length(K))
Qx,Qy,Qz = zeros(length(K)),zeros(length(K)),zeros(length(K))
Ax,Ay,Az = zeros(length(K)),zeros(length(K)),zeros(length(K))

@everywhere i = 0;
for kk in 1:length(K)
    @everywhere i+=1;
    @everywhere k = K[i]
    println("Run $i: Reconstruction gains for k=$(k)")    

    # Define fx ssa function to act on time domian δ for each 
    @everywhere fssa(δ)   = fx_process(δ,dt,fmin,fmax,imp_ssa,  (k))
    @everywhere fqssa(δ)  = fx_process(δ,dt,fmin,fmax,imp_qssa, (k))    
    @everywhere faqssa(δ) = fx_process(δ,dt,fmin,fmax,imp_aqssa,(k))  

    rzx,rzx_diff,xα,rx = pmap_fx_ssa(dx, dzx)
    rzy,rzy_diff,yα,ry = pmap_fx_ssa(dy, dzy)
    rzz,rzz_diff,zα,rz = pmap_fx_ssa(dz, dzz)
   
    (qzx,qzy,qzz),(qzx_diff,qzy_diff,qzz_diff),(qαx,qαy,qαz),(qx,qy,qz) = pmap_fx_qssa(Qt,Qc)
    (azx,azy,azz),(azx_diff,azy_diff,azz_diff),(aαx,aαy,aαz),(ax,ay,az) = pmap_fx_aqssa(Qt,Qc)

    Rx[i] = rx;
    Ry[i] = ry;
    Rz[i] = rz;    

    Qx[i] = qx;
    Qy[i] = qy;
    Qz[i] = qz;    

    Ax[i] = ax;
    Ay[i] = ay;
    Az[i] = az;    

   
    println("Rx=$(rx) Ry=$(ry) and Rz=$(rz)")
    println("Rx=$(qx) Ry=$(qy) and Rz=$(qz)")
    println("Rx=$(ax) Ry=$(ay) and Rz=$(az)")    
    
end