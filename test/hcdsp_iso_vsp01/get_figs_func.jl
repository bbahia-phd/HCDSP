function scalar_get_fig_1(;ls=16,ts=14)
    j = 121; r = 15;
    p=99;

    @show a = max(maximum(vec(ozz)),maximum(vec(dbz)))*0.1

    figname = "time_slice_($j)_receiver_($r)"
    figure(figname,figsize= (10,5))

    subplot(121)
    SeisPlotTX(dzz[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=ts,
               labelsize=ls,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(a)")

    subplot(122)
    SeisPlotTX(dbz[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=ts,
               labelsize=ls,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(b)")

    tight_layout()
#    savefig(joinpath(dev_dir,"files",figname))
end

function get_fig_2()
    j = 121; r = 15;
    p=99;
    a = max(maximum(vec(dzz)),maximum(vec(db)))*0.1

    figname = "res_time_slice_($j)_receiver_($r)"
    figure(figname,figsize= (15,10) .* 1.0)

    subplot(231)
    SeisPlotTX(pgd_fkt[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(a)")

    subplot(232)
    SeisPlotTX(fp_fkt[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(b)")

    subplot(233)
    SeisPlotTX(admm_fkt[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(c)")

    subplot(234)
    SeisPlotTX(r1[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(d)")

    subplot(235)
    SeisPlotTX(r2[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(e)")

    subplot(236)
    SeisPlotTX(r3[j,:,:,r],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(f)")
     
    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end

function single_get_fig_2()
    j = 121; r = 15;
    p=99;
    a = max(maximum(vec(dzz)),maximum(vec(db)))*0.1

    figname = "res_time_slice_($j)_receiver_($r)"
    figure(figname,figsize= (15,10) .* 1.0)

    subplot(231)
    SeisPlotTX(pgd_fkt[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(a)")

    subplot(232)
    SeisPlotTX(fp_fkt[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(b)")

    subplot(233)
    SeisPlotTX(admm_fkt[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(c)")

    subplot(234)
    SeisPlotTX(r1[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(d)")

    subplot(235)
    SeisPlotTX(r2[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(e)")

    subplot(236)
    SeisPlotTX(r3[j,:,:],
               ox = sx_min,
               dx = drz,
               xlabel="Source x  position (km)",
               xticks=[8000,9000,10000],
               xticklabels=["8","9","10"],
               ticksize=14,
               labelsize=16,
               oy = sy_min,
               dy = drz,
               ylabel="Source y  position (km)",
               yticks=[8200,9000,10000,11000],
               yticklabels=["8","9","10","11"],
               cmap="gray",
               vmin=-a,
               vmax=a,
               pclip=p,
               fignum=figname,
               title="(f)")
     
    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end

function get_fig_3(;h=7,w=5)
    sx = 105; sy = 105;
    p=99;
    a = max(maximum(vec(dzz)),maximum(vec(dbz)))*0.1
    xc = 1.8;
    wsc = maximum(vec(abs.(dzz)))*0.2

    figname = "res_csg_($sx)_($sy)"
    figure(figname,figsize= (w,h) .* 1.0)

    subplot(231)
    SeisPlotTX(pgd_fkt[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(a)")

    subplot(232)
    SeisPlotTX(fp_fkt[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(b)")

    subplot(233)
    SeisPlotTX(admm_fkt[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(c)")

    subplot(234)
    SeisPlotTX(r1[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(d)")

    subplot(235)
    SeisPlotTX(r2[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(e)")

    subplot(236)
    SeisPlotTX(r3[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(f)")
     
    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end 

function get_fig_4()
    sx = 105; sy = 105;
    p=99;
    a = max(maximum(vec(dzz)),maximum(vec(db)))*0.1
    xc = 1.8;
    wsc = maximum(vec(abs.(dzz)))*0.2

    figname = "input_csg_($sx)_($sy)"
    figure(figname,figsize= (7,6) .* 1.0)

    subplot(121)
    SeisPlotTX(dzz[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,              
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(a)")

    subplot(122)
    SeisPlotTX(db[:,sx,sy,:],
               ox = -gl_min,
               dx = drz,
               xlabel="Receiver depth (km)",
               xticks=[1400,1600,1800],
               xticklabels=["1.4","1.6","1.8"],
               ticksize=14,
               labelsize=16,
               oy = 0.0,
               dy = dt,
               ylabel="t(s)",
               yticks=[0.0,1.0,2.0],
               yticklabels=["0.0","1.0","2.0"],
               style="wiggles",
               scal = 1/wsc,
               xcur=xc,
#               cmap="gray",
#               vmin=-a,
#               vmax=a,
#               pclip=p,
               fignum=figname,
               title="(b)")
     
    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end 
