pwd()

using Distributed
addprocs(5)

@everywhere dev_dir=joinpath(homedir(),"projects")
@everywhere using Pkg
@everywhere Pkg.activate(joinpath(dev_dir,"HCDSP"))
Pkg.status()

using SeisMain, SeisPlot, PyPlot, StatsBase
using SeisProcessing
using HCDSP

# Set path
work =  joinpath(homedir(),"Desktop/data/seismic_data/valhall/");
bin_files = joinpath(work,"bin");
su_files=joinpath(work,"su");
seis_files=joinpath(work,"seis");
figs=joinpath(work,"figs");

file_name=readdir(su_files)
file_name=file_name[ findall( x -> occursin(".su",x), file_name ) ]

data_path=joinpath.(su_files,file_name[ findall( x -> occursin(".su",x), file_name ) ])
save_path=similar(data_path);

nf=length(data_path);
d=Array{Array{Float32},1}(undef,nf);
h=Array{Array{Header},1}(undef,nf);
ext=Array{Any,1}(undef,nf);

# Read each file
for i in eachindex(data_path)
    indx = findfirst(isequal('.'), file_name[i])-1
#    save_path[i]=joinpath(seis_files,"$(file_name[i][1:indx]).seis")
    save_path[i]="$(file_name[i][1:indx]).seis"
    
    println("Writing $(data_path[i])")
    SegyToSeis(data_path[i],save_path[i],format="su")
    
    dtmp,htmp,etmp=SeisRead(save_path[i])
    d[i]=dtmp;
    h[i]=htmp;
    ext[i]=etmp;    
end

# Set new data, header and extent
nt=2000; n2=0;
for i in 1:nf
    n2 += ext[i].n2
end
din = zeros(Float32,nt,n2);
hin = Vector{Header}(undef,n2);

i1=0;i2=0;
for i in 1:nf
    # indexes
    i1=i2+1;
    i2=i1+ext[i].n2-1
    
    println(size(d[i]))

    # place
    din[:,i1:i2] += d[i]
    hin[i1:i2] .= h[i]
end

# Write seis file (the crg)
extin = ext[1];
extin.n2=n2;
SeisWrite(joinpath(work,"seis/sailline_638-738_july22.seis"),din,hin,extin);

# Read CRG
d,h,ext=SeisRead(joinpath(work,"seis/sailline_638-738_july22.seis"));

sx = SeisMain.ExtractHeader(h,"sx");
sy = SeisMain.ExtractHeader(h,"sy");

gx = SeisMain.ExtractHeader(h,"gx");
gy = SeisMain.ExtractHeader(h,"gy");

x0 = min(sx |> minimum, gx |> minimum)
y0 = min(sy |> minimum, gy |> minimum)

# Set the system's origing
sx .-= x0;
sy .-= y0;
gx .-= x0;
gy .-= y0;

# Estimated before
sr_az = 65;

# cos and sine for rotation
c = cosd(sr_az); s = sind(sr_az);

# rotated source position
sxr = sx .* c .- sy .* s;
syr = sx .* s .+ sy .* c;

# rotated receiver pos
gxr = gx .* c .- gy .* s;
gyr = gx .* s .+ gy .* c;

# find the receivers you want
itr = 0; inds = [];
sxb = []; syb = [];
gxb = []; gyb = [];
for i in 1:ext.n2
    isx = sxr[i]; isy = syr[i];
    igx = gxr[i]; igy = gyr[i];
    if (-9000 < isx < 2e3) && (6000 < isy < 13000)
        itr += 1;
        append!(inds,i);
        append!(sxb,isx);
        append!(syb,isy);
        append!(gxb,igx);
        append!(gyb,igy);
    end
end

n2 = length(inds)
@assert n2 == itr
ext.n2 = n2;   

# for mid points
mx = similar(sxb)
my = similar(syb)

# for offsets
hx = similar(sxb)
hy = similar(syb)

# for azimuth
h  = similar(sxb)
az = similar(syb)

# manually define the SeisHeader for the raw data
F32 = Float32
headjl = Vector{Header}()

# manually extract the traces too
dout = zeros(F32,nt,n2);

# some fields are left blank (0.0)
for n in 1:n2  
    # mid points x and y
    mx[n] = (sxb[n]+gxb[n])/2;
    my[n] = (syb[n]+gyb[n])/2;

    # offsets
    hx[n] = (sxb[n]-gxb[n]);
    hy[n] = (syb[n]-gyb[n]);
    h[n]  = sqrt(hx[n]^2+hy[n]^2);

    az[n] = (180/π) * atan(hy[n],hx[n])
    if az[n] < 0
        az[n] += 360.0
    end
    
    hd = SeisMain.InitSeisHeader()
    hd.tracenum = n
    hd.o1    = 0.f0
    hd.n1    = Int32(nt)
    hd.d1    = F32(0.004)
    hd.sx    = F32(sxb[n])
    hd.sy    = F32(syb[n])
    hd.gx    = F32(gxb[n])
    hd.gy    = F32(gyb[n])
    hd.mx    = F32(mx[n])
    hd.my    = F32(my[n])
    hd.hx    = F32(hx[n])
    hd.hy    = F32(hy[n])
    hd.h     = F32(h[n])
    hd.az    = F32(az[n])

    headjl = push!(headjl,hd)

    # get the right trace
    i = inds[n]
    dout[:,n] .= d[:,i];
end

# define a raw ext file
raw_ext = SeisMain.Extent(nt,n2, 1, 1, 1,               
                          0.f0, 0.f0, 0.f0, 0.f0, 0.f0,
                          0.004f0, 1, 0, 0, 0,
                          "time","sources","","","",
                          "s","index","","","","rotated_valhall")

# Write seis file after band-pass filtering
SeisWrite(joinpath(work,"seis/rotated_valhall_raw.seis"), dout, headjl, raw_ext)

# Read windowed CRG
d,h,e=SeisRead(joinpath(work,"seis/rotated_valhall_raw.seis"));

# clipping the common receiver gather
dc = copy(d);

# Clip to see later reflections a little... NB: This is no good
#val = quantile(abs.(vec(dc)), 0.9999)
#dc = clamp.(dc, -val, val);

# Gain to see later reflections a little... NB: This is also no good.
#dg = SeisGain(dc,dt=0.004,kind="time",coef=[1.1,0.0]); # slow & type unstable

# Band-pass
#f1,f2,f3,f4 = 0,10,45,60;
#dbp = SeisBPFilter(dg,0.004,f1,f2,f3,f4); # slow & type unstable

# F32 
#dc = convert.(Float32,dc);

# Plot
SeisPlotTX(dc[:,1:1000],wbox=10,hbox=3,pclip=90,cmap="gray")

# These are already new coords so I am rewritting stuff
sx = SeisMain.ExtractHeader(h,"sx"); sy = SeisMain.ExtractHeader(h,"sy");
gx = SeisMain.ExtractHeader(h,"gx"); gy = SeisMain.ExtractHeader(h,"gy");

# minimum to set up origin of reg grid
sx_max = maximum(sx);  sx_min = minimum(sx);
sy_max = maximum(sy);  sy_min = minimum(sy);

# grind spacing
dsx = Int32(50); dsy = Int32(50);

# grid size is nsx × nsy
nsx = floor(Int32,(sx_max-sx_min)/dsx)+1
nsy = floor(Int32,(sy_max-sy_min)/dsy)+1

# binned data
dbin = zeros(Float32,nt,nsx,nsy);

# trace counter per bin
count = zeros(Float32,nsx,nsy);

# sampling operator per bin
T = zeros(Float32,nsx,nsy);

# Initialize regular grids
sx_grid = []; sy_grid = [];

for itr in 1:size(dc,2)

    ix=floor(Int32,(sx[itr] - sx_min)./dsx)+1;
    iy=floor(Int32,(sy[itr] - sy_min)./dsy)+1;
  
    dbin[:,ix,iy] .+= dc[:,itr];
    count[ix,iy] += 1.0;
    T[ix,iy]  = elt(1.0);

    tsx = sx_min + (ix-1)*dsx;
    append!(sx_grid, tsx);

    tsy = sy_min + (iy-1)*dsy;
    append!(sy_grid, tsy);
end

# (sx,sy) geometry for CRG
scatter(sxb,syb,label="Grid obs")
scatter(sx_grid,sy_grid,label="Reg Grid")
legend()

# full (desired) grid
sx_full = sx_min:dsx:sx_max |> collect; sx_full = repeat(sx_full,nsy);
sy_full = sy_min:dsy:sy_max |> collect; sy_full = repeat(sy_full,inner=nsx);

# manually define the SeisHeader for the raw data
headjl = Vector{Header}();

for itr in eachindex(T)
    
    ix=floor(Int32,(sx_full[itr] - sx_min)./dsx)+1;
    iy=floor(Int32,(sy_full[itr] - sy_min)./dsy)+1;
  
    tsx = sx_min + (ix-1)*dsx;
    tsy = sy_min + (iy-1)*dsy;

    hd = SeisMain.InitSeisHeader()

    hd.tracenum = itr
    hd.o1    = 0.f0
    hd.n1    = Int32(nt)
    hd.d1    = F32(0.004)
    hd.sx    = F32(tsx)
    hd.sy    = F32(tsy)
    hd.gx    = F32(gxb[1])
    hd.gy    = F32(gyb[1])

    headjl = push!(headjl,hd)
end

# define a raw ext file
raw_ext = SeisMain.Extent(nt,nsx, nsy, 1, 1,               
                          0.f0, sx_min, sy_min, 0.f0, 0.f0,
                          0.004f0, dsx, dsy, 0.f0, 0.f0,
                          "time","sx","sy","","",
                          "s","m","m","","","binned_rotated_valhall")

# Write seis file after binning
SeisWrite(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_$(dsx)x$(dsy).seis"), dbin, headjl, raw_ext)
#SeisWrite(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_50x50.seis"), dbin, headjl, raw_ext)

# Read seis file after binning
di,hi,ei = SeisRead(joinpath(work,"seis/binned_rotated_valhall_sailline_638-738_$(dsx)x$(dsy).seis"));

# binned data size
@assert dbin == di

# average repeated bins
cindex = CartesianIndices((1:nsx,1:nsy));
tot = 0; tot_bin = prod((nsx,nsy));
for i in eachindex(T)
    if T[i] == 1
        tot += 1;
        dbin[:,cindex[i]] ./= count[i]
    end
end

# percentage of alive traces
trc_perc = round(tot/tot_bin*100,digits=2)
println("$tot out of $tot_bin ($trc_perc %) alive traces")

# percentage for fig clipping
pc=90;

# inline
j = 50;
SeisPlotTX(dbin[:,j,:],pclip=pc,cmap="gray")

# xline
j = 30;
SeisPlotTX(dbin[:,:,j],pclip=pc,cmap="gray")

# Write binaries
read_write(joinpath(work,"bin/binned_rotated_sailline_638-738_50x50.bin"),
           "w" ; input=dbin, T = elt );
read_write(joinpath(work,"bin/binned_sampling_sailline_638-738_50x50.bin"),"w"
           ; input=T, T = elt );

#read_write(joinpath(work,"bin/binned_rotated_sailline_638-738_100x100.bin"),  "w" ; input=dbin );
#read_write(joinpath(work,"bin/sampling_sailline_638-738_100x100.bin"), "w" ; input=T );
