function pmap_fx_ssa(δ,δc)
    # apply patching on input
    in_patch,pid   = fwdPatchOp(δ, psize, polap, smin, smax);

    # SSA all patches
    out_patch  = similar(in_patch);
    out_patch .= pmap(fssa,in_patch);    

    # un-patch
    out_full = adjPatchOp(out_patch, pid, psize, polap, smin, smax);

    # scale out to clean
    α = dot(δc,out_full) / dot(out_full,out_full)
    
    # get differences
    diff_full  = δc .- α .* out_full;

    # quality
    r = quality(δc,out_full)

    return out_full,diff_full,α,r
end

function pmap_fx_qssa(δ,δc)

    # componentwise inputs
    xin = imagi.(δc)
    yin = imagj.(δc)
    zin = imagk.(δc)

    # apply patching on input
    in_patch, pid   = fwdPatchOp(δ, psize, polap, smin, smax);

    # SSA all patches
    out  = similar(in_patch);
    out .= pmap(fqssa,in_patch);    

    # un-patch
    xout = adjPatchOp(imagi.(out), pid, psize, polap, smin, smax);
    yout = adjPatchOp(imagj.(out), pid, psize, polap, smin, smax);
    zout = adjPatchOp(imagk.(out), pid, psize, polap, smin, smax);

    # scale out to clean
    αx = dot(xin,xout) / dot(xout,xout)
    αy = dot(yin,yout) / dot(xout,xout)
    αz = dot(zin,zout) / dot(xout,xout)
    
    # get differences
    diffx  = xin .- αx .* xout;    
    diffy  = yin .- αy .* yout;    
    diffz  = zin .- αz .* zout;    

    # quality
    rx = quality(xin,xout)
    ry = quality(yin,yout)
    rz = quality(zin,zout)

    return (xout,yout,zout),(diffx,diffy,diffz),(αx,αy,αz),(rx,ry,rz)
end

function pmap_fx_aqssa(δ,δc)

    # componentwise inputs
    xin = imagi.(δc)
    yin = imagj.(δc)
    zin = imagk.(δc)

    # apply patching on input
    in_patch, pid   = fwdPatchOp(δ, psize, polap, smin, smax);

    # SSA all patches
    out  = similar(in_patch);
    out .= pmap(faqssa,in_patch);    

    # un-patch
    xout = adjPatchOp(imagi.(out), pid, psize, polap, smin, smax);
    yout = adjPatchOp(imagj.(out), pid, psize, polap, smin, smax);
    zout = adjPatchOp(imagk.(out), pid, psize, polap, smin, smax);

    # scale out to clean
    αx = dot(xin,xout) / dot(xout,xout)
    αy = dot(yin,yout) / dot(xout,xout)
    αz = dot(zin,zout) / dot(xout,xout)
    
    # get differences
    diffx  = xin .- αx .* xout;    
    diffy  = yin .- αy .* yout;    
    diffz  = zin .- αz .* zout;    

    # quality
    rx = quality(xin,xout)
    ry = quality(yin,yout)
    rz = quality(zin,zout)

    return (xout,yout,zout),(diffx,diffy,diffz),(αx,αy,αz),(rx,ry,rz)
end
