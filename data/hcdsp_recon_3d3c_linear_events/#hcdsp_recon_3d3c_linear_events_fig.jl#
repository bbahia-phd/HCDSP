pwd()

cd("./hcdsp_recon_3d3c_linear_events")

using Pkg
Pkg.activate("./../../")
Pkg.status()

using Revise
using HDF5
using PyPlot
using SeisMain, SeisPlot
using StatsBase,Statistics
using HCDSP

nt=100;n1=40;n2=40;


bin_files = filter(x->occursin(".bin",x),readdir())

function reread(fname)
    global nt,n1,n2
    out = read_write(fname,"r",n=(nt,n1,n2),T=Float32)
    out = reshape(out,(nt,n1,n2));
    return out
end

# ideal
dzx = reread(bin_files[13]);
dzy = reread(bin_files[14]);
dzz = reread(bin_files[15]);

# noisy
dnx = reread(bin_files[4]);
dny = reread(bin_files[5]);
dnz = reread(bin_files[6]);

# ssa
rzx = reread(bin_files[10]);
rzy = reread(bin_files[11]);
rzz = reread(bin_files[12]);

# qssa
qzx = reread(bin_files[7]);
qzy = reread(bin_files[8]);
qzz = reread(bin_files[9]);

# aqssa
azx = reread(bin_files[1]);
azy = reread(bin_files[2]);
azz = reread(bin_files[3]);

function getfigs()

   close("all")

    fig1()

    fig2()

end

function fig1(;cmap="gray",y=10,xc=2)
    # this means 10*10 = 100m offset in x
    fig_num="offsety_100m"

    figure(fig_num,figsize=(6,6.5))

    subplot(231);
    SeisPlotTX(dzx[:,1:2:end,y],fignum=fig_num,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(a)")

    subplot(232);
    SeisPlotTX(dzy[:,1:2:end,y],fignum=fig_num,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(b)")

    subplot(233);
    SeisPlotTX(dzz[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(c)")

    subplot(234);
    SeisPlotTX(dnx[:,1:2:end,y],fignum=fig_num, 
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(d)")

    subplot(235);
    SeisPlotTX(dny[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(e)")

    subplot(236);
    SeisPlotTX(dnz[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(f)")

    tight_layout()
end
 
function fig2(;cmap="gray",y=10)
    
    fig_num="offsety_100m_results"

    figure(fig_num,figsize=(6,12))

    subplot(331);
    SeisPlotTX(rzx[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(a)")
               

    subplot(332);
    SeisPlotTX(rzy[:,1:2:end,y],fignum=fig_num,
              xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(b)")
               
    subplot(333);
    SeisPlotTX(rzz[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(c)")
               

    subplot(334);
    SeisPlotTX(qzx[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(d)")
                          

    subplot(335);
    SeisPlotTX(qzy[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(e)")
               
    subplot(336);
    SeisPlotTX(qzz[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(f)")
               
               
    subplot(337);
    SeisPlotTX(azx[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(g)")
               
               
    subplot(338);
    SeisPlotTX(azy[:,1:2:end,y],fignum=fig_num,
              xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(h)")
               

    subplot(339);
    SeisPlotTX(azz[:,1:2:end,y],fignum=fig_num,
               xcur=2,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(i)")              
               
    tight_layout()
    
end

function geta(ideal,approx)
    return dot(vec(ideal),vec(approx))/dot(vec(approx),vec(approx))
end

function fig3(;cmap="gray",y=10,xc=1,sc=5)
    
    fig_num="offsety_100m_results_difference"

    figure(fig_num,figsize=(6,12))

    ax = maximum(vec(dzx));
    bx = -ax;

    ay = maximum(vec(dzy));
    by = -ay;

    az = maximum(vec(dzz));
    bz = -az;

    αx = geta(dzx,rzx);
    αy = geta(dzy,rzy);
    αz = geta(dzz,rzz);
    
    dx = sc .* (dzx .- αx .* rzx);
    dy = sc .* (dzy .- αy .* rzy);
    dz = sc .* (dzz .- αz .* rzz);
    
    subplot(331);
    SeisPlotTX(dx[:,1:2:end,y],fignum=fig_num,
               scal=1/ax,
               vmin=bx,
               vmax=ax,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(a)")
               
    subplot(332);
    SeisPlotTX(dy[:,1:2:end,y],fignum=fig_num,
               scal=1/ay,
               vmin=by,
               vmax=ay,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(b)")
               
    subplot(333);
    SeisPlotTX(dz[:,1:2:end,y],fignum=fig_num,
               scal=1/az,
               vmin=bz,
               vmax=az,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(c)")
               

    αx = geta(dzx,qzx);
    αy = geta(dzy,qzy);
    αz = geta(dzz,qzz);

    dx = sc .* (dzx .- αx .* qzx);
    dy = sc .* (dzy .- αy .* qzy);
    dz = sc .* (dzz .- αz .* qzz);

    subplot(334);
    SeisPlotTX(dx[:,1:2:end,y],fignum=fig_num,
               scal=1/ax,
               vmin=bx,
               vmax=ax,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(d)")
                          

    subplot(335);
    SeisPlotTX(dy[:,1:2:end,y],fignum=fig_num,
               scal=1/ay,
               vmin=by,
               vmax=ay,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(e)")
               
    subplot(336);
    SeisPlotTX(dz[:,1:2:end,y],fignum=fig_num,
               scal=1/az,
               vmin=bz,
               vmax=az,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(f)")

    αx = geta(dzx,azx);
    αy = geta(dzy,azy);
    αz = geta(dzz,azz);

    dx = sc .* (dzx .- αx .* azx);
    dy = sc .* (dzy .- αy .* azy);
    dz = sc .* (dzz .- αz .* azz);
               
    subplot(337);
    SeisPlotTX(dx[:,1:2:end,y],fignum=fig_num,
               scal=1/ax,
               vmin=bx,
               vmax=ax,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(g)")
               
               
    subplot(338);
    SeisPlotTX(dy[:,1:2:end,y],fignum=fig_num,
               scal=1/ay,
               vmin=by,
               vmax=ay,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(h)")
               

    subplot(339);
    SeisPlotTX(dz[:,1:2:end,y],fignum=fig_num,
               scal=1/az,
               vmin=bz,
               vmax=az,
               xcur=xc,
               dy=0.004,
               yticks=[0,0.1,0.2,0.3],
               yticklabels=[0,0.1,0.2,0.3],
               ylabel="Time (s)",
               xlabel="Offset x (m)",
               dx=20,
               style="overlay",
               cmap=cmap,
               title="(i)")              
               
    tight_layout()
    
end
