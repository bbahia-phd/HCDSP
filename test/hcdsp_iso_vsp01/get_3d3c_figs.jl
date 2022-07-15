function vector_get_fig_1(;ls=18,ts=16)
    j = 91; r = 15;
    p=99;

    @show a = max(maximum(vec(dbx)),maximum(vec(dby)),maximum(vec(dbz)))*0.1

    figname = "3c_time_slice_($j)_receiver_($r)"
    figure(figname,figsize=(15,10))

    subplot(231)
    SeisPlotTX(dzx[j,:,:],
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
    SeisPlotTX(dzy[j,:,:],
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
    SeisPlotTX(dzz[j,:,:],
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
    SeisPlotTX(dbx[j,:,:],
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
    SeisPlotTX(dby[j,:,:],
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
    SeisPlotTX(dbz[j,:,:],
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

function vector_get_fig_2(;w=15,h=10,ts=16,ls=18)
    j = 91; r = 15;
    p=99;

    figname = "qssa_all_time_slice_($j)_receiver_($r)"

    @show a = max(maximum(vec(dbx)),maximum(vec(dby)),maximum(vec(dbz)))*0.1

    figure(figname,figsize= (w,h) .* 1.0)

    subplot(331)
    SeisPlotTX(rzx[j,:,:],
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
    SeisPlotTX(rzy[j,:,:],
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
    SeisPlotTX(rzz[j,:,:],
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
    SeisPlotTX(qzx[j,:,:],
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
    SeisPlotTX(qzy[j,:,:],
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
    SeisPlotTX(qzz[j,:,:],
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
    SeisPlotTX(azx[j,:,:],
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
    SeisPlotTX(azy[j,:,:],
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
    SeisPlotTX(azz[j,:,:],
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


function vector_get_fig_3(;w=15,h=10,ts=16,ls=18)
    j = 121; r = 15;
    p=99;

    figname = "qssa_all_res_time_slice_($j)_receiver_($r)"

    @show a = max(maximum(vec(dbx)),maximum(vec(dby)),maximum(vec(dbz)))*0.5

    figure(figname,figsize= (w,h) .* 1.0)

   
    @show α = dot(dzx,rzx) / dot(rzx,rzx);
    p1 = dzx[j,:,:] .- rzx[j,:,:];
    @show α = dot(dzy,rzy) / dot(rzy,rzy);
    p2 = dzy[j,:,:] .- rzy[j,:,:];
    @show α = dot(dzz,rzz) / dot(rzz,rzz);
    p3 = dzz[j,:,:] .- α .* rzz[j,:,:];

    @show α = dot(dzx,qzx) / dot(qzx,qzx);
    p4 = dzx[j,:,:] .- α .* qzx[j,:,:];
    @show α = dot(dzy,qzy) / dot(qzy,qzy);
    p5 = dzy[j,:,:] .- α .* qzy[j,:,:];
    @show α = dot(dzz,qzz) / dot(qzz,qzz);
    p6 = dzz[j,:,:] .- α .* qzz[j,:,:];

    @show α = dot(dzx,azx) / dot(azx,azx);
    p7 = dzx[j,:,:] .- α .* azx[j,:,:];
    @show α = dot(dzy,azy) / dot(azy,azy);
    p8 = dzy[j,:,:] .- α .* azy[j,:,:];
    @show α = dot(dzz,azz) / dot(azz,azz);
    p9 = dzz[j,:,:] .- α .* azz[j,:,:];

    
    subplot(331)
    SeisPlotTX(p1,
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
               title="(a) "*L"{\bf U}_{x} " *" (SSA)")

    subplot(332)
    SeisPlotTX(p2,
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
               title="(b) "*L"{\bf U}_{y} " *" (SSA)")

    subplot(333)
    SeisPlotTX(p3,
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
               title="(c) "*L"{\bf U}_{z} " *" (SSA)")

    subplot(334)
    SeisPlotTX(p4,
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
               title="(d) "*L"{\bf U}_{x} " *" (QSSA)")


    subplot(335)
    SeisPlotTX(p5,
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
               title="(e) "*L"{\bf U}_{y} " *" (QSSA)")

    subplot(336)
    SeisPlotTX(p6,
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
               title="(f) "*L"{\bf U}_{z} " *" (QSSA)")
     
    subplot(337)
    SeisPlotTX(p7,
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
                title="(g) "*L"{\bf U}_{x} " *" (AQSSA)")

    subplot(338)
    SeisPlotTX(p8,
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
                title="(h) "*L"{\bf U}_{y} " *" (AQSSA)")

    subplot(339)
    SeisPlotTX(p9,
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
                title="(i) "*L"{\bf U}_{z} " *" (AQSSA)")

    tight_layout()
    savefig(joinpath(dev_dir,"files",figname))
end
