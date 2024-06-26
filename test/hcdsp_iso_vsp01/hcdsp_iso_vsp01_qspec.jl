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
data_path = "/media/bbahia/DATA/seismic_data/iso_vsp01/crg1350";  

dzx,hzx,ext_zx = SeisRead(joinpath(data_path,"iso_vsp01_zx_crg1350.seis"));
dzy,hzy,ext_zy = SeisRead(joinpath(data_path,"iso_vsp01_zy_crg1350.seis"));
dzz,hzz,ext_zz = SeisRead(joinpath(data_path,"iso_vsp01_zz_crg1350.seis"));

dt = Float64(ext_zz.d1);
n = nt,n1,n2 = 217,205,205;
drz = 16.6;

dzx = permutedims(reshape(dzx,n),(1,3,2));
dzy = permutedims(reshape(dzy,n),(1,3,2));
dzz = permutedims(reshape(dzz,n),(1,3,2));

rz_init = 1350.0; # m
rz_end  = 1850.01; # m
dsx = dsy = drz = 16.6 ; # m 
rz_axis = range(rz_init, rz_end, step=drz);
nr = length(rz_axis)

# get geometry
sx = SeisMain.ExtractHeader(hzx,"sx");
sy = SeisMain.ExtractHeader(hzx,"sy");
gl = SeisMain.ExtractHeader(hzx,"gelev");

sx_n,sy_n = sx |> copy,sy |> copy;

# grid
sx_min = sx |> minimum;
sx_max = sx |> maximum;

sy_min = sy |> minimum;
sy_max = sy |> maximum;

gl_min = gl |> maximum;
gl_max = gl |> minimum;

j = 50;

dx = Float64.(dzx[:,:,j]);
dy = Float64.(dzy[:,:,j]);
dz = Float64.(dzz[:,:,j]);

cmap_time = "gray"
cmap_fk   = "jet"
fignum    = "3c_time_domain";
ylabel    = "t"
yunit     = "(s)"
xlabel    = "Source "*L"x"*" position"
xunit     = "(km)"
ts=14;
ls=16;

clf(); close("all")
# Time domain 2D plot #
figure(fignum,figsize=(12,10))

subplot(2,3,1)
SeisPlotTX(dx,
        ox = sx_min,
        dx = drz,
        xticks=[8000,9000,10000],
        xticklabels=["8","9","10"],
        ticksize=ts,
        labelsize=ls,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(a)")

subplot(2,3,2)
SeisPlotTX(dy,
        ox = sx_min,
        dx = drz,
        xticks=[8000,9000,10000],
        xticklabels=["8","9","10"],
        ticksize=ts,
        labelsize=ls,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(b)")

subplot(2,3,3)
SeisPlotTX(dz,
        ox = sx_min,
        dx = drz,
        xticks=[8000,9000,10000],
        xticklabels=["8","9","10"],
        ticksize=ts,
        labelsize=ls,
        oy=0.0,
        dy=dt,
        fignum=fignum,
        cmap=cmap_time,
        ylabel=ylabel,
        yunits=yunit,
        xlabel=xlabel,
        xunits=xunit,
        title="(c)")

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
        title="(d)",
        ticksize=ts,
        labelsize=ls
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
        title="(e)",
        ticksize=ts,
        labelsize=ls
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
        title="(f)",
        ticksize=ts,
        labelsize=ls
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
