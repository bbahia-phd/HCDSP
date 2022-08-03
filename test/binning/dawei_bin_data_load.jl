import Pkg
Pkg.activate(joinpath(homedir(),"projects/HCDSP"))
Pkg.status()

using Revise
using PyPlot

using MAT, MATLAB
mat"addpath('./../../../QSEIS/src/Plotting/')"

using SeisMain, SeisPlot, SeisProcessing

ENV["DATAPATH"]=""

file = matopen("data_patch.mat")
varnames = keys(file)
head = read(file,"H");

dt  = head["dt"][1];
cdp = head["cdp"][1,:]  |> Vector{Float64};
sx  = head["sx"][1,:]   |> Vector{Float64};
sy  = head["sy"][1,:]   |> Vector{Float64};

gx = head["gx"][1,:]   |> Vector{Float64};
gy = head["gy"][1,:]   |> Vector{Float64};

dout = read(file,"D");
dout = SeisBPFilter(dout,dt / 1e6 ,1,3,50,65);

nt,ntrace = size(dout);
t = collect((0:nt-1) .* dt) ./ 1e6;

x0 = 0.0#min(minimum(sx),minimum(gx)) + 10
y0 = 0.0#min(minimum(sy),minimum(gy)) + 10

sx .-= x0
sy .-= y0
gx .-= x0
gy .-= y0

sr_az = 39π/180;
c = cos(sr_az); s = sin(sr_az);

sxr = sx .* c .- sy .* s;
syr = sx .* s .+ sy .* c;

gxr = gx .* c .- gy .* s;
gyr = gx .* s .+ gy .* c;

#scatter(gx,gy)
#scatter(gxr,gyr)

# for mid points
mx = similar(sxr)
my = similar(syr)

# for offsets
hx = similar(sxr)
hy = similar(syr)

# for azimuth
h = similar(sxr)
az = similar(syr)

# manually define the SeisHeader for the raw data
F32 = Float32
headjl = Vector{Header}()

# some fields are left blank (0.0)
for n in 1:ntrace
   
    # mid points x and y
    mx[n] = (sx[n]+gx[n])/2;
    my[n] = (sy[n]+gy[n])/2;

    # offsets
    hx[n] = (sx[n]-gx[n]);
    hy[n] = (sy[n]-gy[n]);
    h[n]  = sqrt(hx[n]^2+hy[n]^2);

    az[n] = (180/π) * atan(hy[n],hx[n])
    if az[n] < 0
        az[n] += 360.0
    end
    
    hd = SeisMain.InitSeisHeader()
    hd.tracenum = n
    hd.o1    = 0.f0
    hd.n1    = Int32(head["ns"][n])
    hd.d1    = F32(head["dt"][n])
    hd.sx    = F32(head["sx"][n])
    hd.sy    = F32(head["sy"][n])
    hd.gx    = F32(head["gx"][n])
    hd.gy    = F32(head["gy"][n])
    hd.mx    = F32(mx[n])
    hd.my    = F32(my[n])
    hd.hx    = F32(hx[n])
    hd.hy    = F32(hy[n])
    hd.h     = F32(h[n])
    hd.az    = F32(az[n])

    headjl = push!(headjl,hd)
end

# define a raw ext file
raw_ext = SeisMain.Extent(nt,ntrace, 1, 1, 1,               
                          0.f0, 0.f0, 0.f0, 0.f0, 0.f0,
                          dt, 1, 0, 0, 0,
                          "time","mx","my","offset","azimuth",
                          "s","m","m","m","deg","raw_data")

SeisWrite("data_bp_raw", dout, headjl, raw_ext)
SeisWindow("data_bp_raw", "data_bp_raw_1s", key=["t"], minval=[799*dt], maxval=[1399*dt])

# Define a regular grid for binning
mx_min = minimum(mx); mx_max = maximum(mx);
my_min = minimum(my); my_max = maximum(my);
h_max = maximum(h); h_min = minimum(h);
az_min = 0; az_max =  2*180;

dmx = 200;
nx = floor(Int,(mx_max-mx_min)/dmx) + 1;

dmy = 200;
ny = floor(Int,(my_max-my_min)/dmy) + 1;

dh = 400;
nh = floor(Int,(h_max - h_min)/dh) + 1;

daz = 45;
naz = floor(Int,(az_max-az_min)/daz)+1;

param1 = Dict(:osx=>  575103, :osy=>5.84003e6,
              :ogx=>  574774, :ogy=>5.847986e6,
              :dmx => dmx,    :dmy => dmy,
              :dh  => dh,     :daz => daz)

SeisGeometry("data_bp_raw_1s"; param1...)

db,hb,eb = SeisRead("data_bp_raw_1s");
