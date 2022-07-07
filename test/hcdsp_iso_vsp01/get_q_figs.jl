function vector_get_fig_1(;ls=18,ts=16)
    j = 121; r = 15;
    p=99;

    @show a = max(maximum(vec(dbx)),maximum(vec(dby)),maximum(vec(dbz)))*0.1

    figname = "3c_time_slice_($j)_receiver_($r)"
    figure(figname,figsize=(15,10))

    subplot(231)
    SeisPlotTX(ozx[j,:,:,r],
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

    subplot(232)
    SeisPlotTX(ozy[j,:,:,r],
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

    subplot(233)
    SeisPlotTX(ozz[j,:,:,r],
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
                title="(c)")

    subplot(234)
    SeisPlotTX(dbx[j,:,:,r],
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
                title="(d)")

    subplot(235)
    SeisPlotTX(dby[j,:,:,r],
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
                title="(e)")

    subplot(236)
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
                title="(f)")                

    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end

function vector_get_fig_2()
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

function vector_get_fig_3(;w=15,h=10,ts=16,ls=18)
    j = 121; r = 15;
    p=99;

    figname = "qssa_all_res_time_slice_($j)_receiver_($r)"

    figure(figname,figsize= (w,h) .* 1.0)
    a = maximum(vec(dzz))*0.2

    subplot(331)
    SeisPlotTX(x_pgd_fkt[j,:,:,r],
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
               title="(a) "*L"{\bf U}_{x} " *" (PGD QFKT)")

    subplot(332)
    SeisPlotTX(y_pgd_fkt[j,:,:,r],
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
               title="(b) "*L"{\bf U}_{y} " *" (PGD QFKT)")

    subplot(333)
    SeisPlotTX(z_pgd_fkt[j,:,:,r],
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
               title="(c) "*L"{\bf U}_{z} " *" (PGD QFKT)")

    subplot(334)
    SeisPlotTX(x_fp_fkt[j,:,:,r],
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
               title="(d) "*L"{\bf U}_{x} " *" (RED-FP QFKT)")


    subplot(335)
    SeisPlotTX(y_fp_fkt[j,:,:,r],
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
               title="(e) "*L"{\bf U}_{y} " *" (RED-FP QFKT)")

    subplot(336)
    SeisPlotTX(z_fp_fkt[j,:,:,r],
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
               title="(f) "*L"{\bf U}_{z} " *" (RED-FP QFKT)")
     
    subplot(337)
    SeisPlotTX(x_admm_fkt[j,:,:,r],
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
                title="(g) "*L"{\bf U}_{x} " *" (RED-ADMM QFKT)")

    subplot(338)
    SeisPlotTX(y_admm_fkt[j,:,:,r],
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
                title="(h) "*L"{\bf U}_{y} " *" (RED-ADMM QFKT)")

    subplot(339)
    SeisPlotTX(z_admm_fkt[j,:,:,r],
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
                title="(i) "*L"{\bf U}_{z} " *" (RED-ADMM QFKT)")

    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end


function vector_get_fig_4(;h=7,w=5,ls=18,ts=16)
    p=99;
    sx = 105; sy = 105;
    xc = 1.8;
    wsc = maximum(vec(abs.(dzz)))*0.2

    figname = "3c_res_csg_($sx)_($sy)"
    figure(figname,figsize= (w,h) .* 1.0)

    subplot(231)
    SeisPlotTX(z_pgd_fkt[:,sx,sy,:],
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
    SeisPlotTX(z_fp_fkt[:,sx,sy,:],
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
    SeisPlotTX(z_admm_fkt[:,sx,sy,:],
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
