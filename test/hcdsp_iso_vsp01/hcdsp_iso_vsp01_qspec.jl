
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

# data dir home
data_path = "./files/vsp3d9c/iso_vsp01/";  

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz_crg1350.seis"));

dt = Float64(ext_zz.d1);
nt = Int64(ext_zz.n1);

ntr = ext_zz.n2;

rz_init = 1350.0; # m
rz_end  = 1850.01; # m
drz = 16.667 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr = length(rz_axis)

nsline = 205; # number of source lines
ns = 205; # number of sources within lines

@assert ns*nsline == ntr == size(dzz,2)

dzx = reshape(dzx,(nt,ns,nsline));
dzy = reshape(dzy,(nt,ns,nsline));
dzz = reshape(dzz,(nt,ns,nsline));

j = 50;

dx = Float64.(dzx[:,:,j]);
dy = Float64.(dzy[:,:,j]);
dz = Float64.(dzz[:,:,j]);

# Add noise
dnx = SeisAddNoise(dx, -2.0, db=true, L=3);
dny = SeisAddNoise(dy, -2.0, db=true, L=3);
dnz = SeisAddNoise(dz, -2.0, db=true, L=3);

cmap_time = "gray"
cmap_fk   = "jet"
fignum    = "3c_time_domain";
ylabel    = "Time"
yunit     = "(s)"
xlabel    = "Trace"
xunit     = ""

clf(); close("all")

# Time domain 2D plot #
figure(fignum,figsize=(10,8))

subplot(2,3,1)
SeisPlotTX(dx,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(a)"
)

subplot(2,3,2)
SeisPlotTX(dy,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(b)"
)

subplot(2,3,3)
SeisPlotTX(dz,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(c)"
)

subplot(2,3,4)
SeisPlotFK(dx,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_fk,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(d)"
)
colorbar()

subplot(2,3,5)
SeisPlotFK(dy,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_fk,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(e)"
)
colorbar()

subplot(2,3,6)
SeisPlotFK(dz,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_fk,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(f)"
)
colorbar()

tight_layout()
gcf()

savefig(fignum,dip=300)

## New fig: Q-SPECTRA ##

fignum="qft_spectra"

# Temporary Quaternion
Qt = quaternion(dx,dy,dz);

μ = μ0

# Q-Amplitude Spectra #

Qft = qfft(Qt,μ,"left",1)

qas = fftshift(sum(abs.(Qft),dims=2))

# Q-FK Spectra #
cmap="PuOr"
pclip=99.9

aspect="auto"
interpolation="Hanning"

titlesize=16
labelsize=14

dpi=300

xlabel = "Wavenumber"
xunits = "(1/m)"

ylabel = "Frequency"
yunits = "Hz"

dk = 1/drz/size(Qt,2)
kmin = -dk*size(Qt,2)/2
kmax =  dk*size(Qt,2)/2

df = 1/dt/size(Qt,1)
fmin = -df*size(Qt,1)/2
fmax =  df*size(Qt,1)/2
nf = Int32(floor(size(Qt,1)/2))
faxis = range(fmin,fmax;length=2nf+1)

Qf = qfft(Qt,μ,"left")

qfk = abs.(fftshift(Qf));

clf(); close("all")
# Time domain 2D plot #
fignum="qft_spectra"
figure(fignum,figsize=(8,4))

subplot(1,2,1)
plot(faxis,qas)
grid()
PyPlot.title("(a)")
PyPlot.xlabel("Frequency (Hz)")
PyPlot.ylabel("Amplitude")

subplot(1,2,2)

imshow(qfk,
        cmap=cmap_fk,
        extent=[kmin,kmax,fmin,fmax],
        aspect=aspect,
        interpolation=interpolation)       
colorbar()
PyPlot.title("(b)")
PyPlot.xlabel("Wavenumber (1/m)")
PyPlot.ylabel("Frequency Hz")

tight_layout()

gcf()

savefig(fignum,dpi=300)
