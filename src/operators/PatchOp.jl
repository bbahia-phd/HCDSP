"""
    patches,pind = fwdPatchOp(in, psize, polap, smin, smax)

Divide `in` in overlapping patches of size `psize` and overlap `olap`.

Returns `patches` and `pind`. `patches` is a `Vector{typeof(in)}` that carries the overlapping patches.
`pind` are the indexes for each patch which is used by `adjPatchOp` for merging the paches.

# Arguments
- `in::AbstractArray{T,N}`: input seismic data as ND array
- `psize::NTuple{N,Int}`: N-tuple containing the patch size
- `polap::NTuple{N,Int}`: N-tuple containing the patch overlap in %
- `smin::NTuple{N,Int}`: N-tuple containing the first index to be patched
- `smax::NTuple{N,Int}`: N-tuple containing the last index to be patched 

# Examples
```julia-repl
julia> in = randn(Float64,nt,nx,ny)
julia> psize, polap, smin, smax = (pt,px,py), (perct,percx,percy), (ot,ox,oy), (ft,fx,fy)
julia> patches,pid = fwdPatchOp(in, psize, polap, smin, smax)
julia> size(patches)
```
""" 
function fwdPatchOp(IN::AbstractArray{T,N},
                    psize::NTuple{N,Int},
                    polap::NTuple{N,Int},
                    smin::NTuple{N,Int},
                    smax::NTuple{N,Int}) where {T <: Number,N}
    
    # case flag == "fwd"

    # actual input dimensions
    in_dim = smax .- smin .+ 1;

    # convert overlaps to samples
    solap = round.(Int64,psize.*polap./100)

    # safe-guard for window bigger than data
    psize = min.(psize , in_dim)

    # number of windows per dimension (nwd) and total (nwt)
    nwd = floor.(Int64,in_dim ./ (psize .- solap))
    nwd = nwd .+ Int.( (smin .+ (nwd .-1).*(psize.-solap) .+ psize) .< smax )
    nwt = prod(nwd)

    # allocate output
    OUT = Array{ Array{T,N}, 1 }(undef,nwt);
    pind = Array{ NTuple{ N, UnitRange{Int} } }(undef,nwt)
    cind = CartesianIndices(nwd);

    for i in eachindex(OUT)
        mini = smin .+ (Tuple(cind[i]) .- 1) .* (psize .- solap)
        maxi = min.(mini .+ (psize .- 1), smax)
        pind[i] = map((x,y) -> x:y, mini, maxi)
        @inbounds OUT[i] = IN[CartesianIndices((pind[i]))]
    end

    return OUT,pind
    
end
#===========================================================================================#

"""
    out = adjPatchOp(patches, pind, psize, polap, smin, smax)

Returns full (merged) section from `patches` and `pind` indexes.

# Arguments
- `patches::AbstractArray{Array{T},N}`: array containing patches of original seismic data
- `pind::Vector{NTuple{N,UnitRange{Int}}}`: corresponding indexes for patching location
- `psize::NTuple{N,Int}`: N-tuple containing the patch size
- `polap::NTuple{N,Int}`: N-tuple containing the patch overlap in %
- `smin::NTuple{N,Int}`: N-tuple containing the first index to be patched
- `smax::NTuple{N,Int}`: N-tuple containing the last index to be patched 

# Examples
```julia-repl
julia> patches,pid = fwdPatch(in, psize, polap, smin, smax)
julia> out = fwdPatchOp(patches, pid, psize, polap, smin, smax)
julia> size(out) == size(in)
```
""" 
function adjPatchOp(IN::AbstractArray{Array{T,N}},
                    pind::Vector{NTuple{N,UnitRange{Int}}},
                    psize::NTuple{N,Int},
                    polap::NTuple{N,Int},
                    smin::NTuple{N,Int},
                    smax::NTuple{N,Int}) where {T <: Number, N}

    # case: flag == "adj"

    #total number of windows (nwt)
    nwt = length(IN)

    # output dimensions
    outd = smax .- smin .+ 1

    # convert overlaps to samples
    solap = round.(Int,psize.*polap./100)

    # allocate output
    OUT = zeros(T,outd)

    for i in eachindex(IN)
        
        tmp = copy(IN[i])
        
        # tapers
        tapi = solap .* Int.( map(x -> x[1]  , pind[i]) .> smin )
        tapf = solap .* Int.( map(x -> x[end], pind[i]) .< smax )
        
        # apply taper
        applyTaper!(tmp,tapi,tapf)
        #applyTaper!(IN[i],tapi,tapf)

        OUT[CartesianIndices((pind[i]))] += tmp;
    end

    return OUT

end
#===========================================================================================#

"""
    applyTaper!(IN,tapi,tapf)

Apply tapers `tapi` and `tapf` to `IN`.
"""
function applyTaper!(IN::Array{T,N},tapi::NTuple{N,Int}, tapf::NTuple{N,Int}) where {T <: Number, N}

    # set of indexes
    cwin = CartesianIndices(IN);
       
    # patches might vary in size
    szt = size(IN);

    for j in cwin
        tap = getTaper.(Tuple(j),szt,tapi,tapf)
        IN[j] *= T.(prod(tap))
    end
end
#===========================================================================================#

"""
    getTaper(ind, n, tapi, tapf)

Obtains tapering values for given index `ind`.
"""
function getTaper(ind::Int, n::Int, tapi::Int, tapf::Int)

    t = 1.0;
    
    if ind >= 1 && ind <= tapi ;
        t = t / (tapi-1) * (ind-1);
    elseif ind > tapi && ind <= n-tapf
        t = t
    elseif ind > n-tapf && ind <= n
        t = t - t / (tapf-1) * (ind-1-n+tapf)
    end

    return t
 
end
#===========================================================================================#

#=
# TODO: Follow convention SeisPatchOp() 

# Multiple dispatch
SeisPatchOp(in,psize,polap,smin,smax)      = fwdPatchOp(in,psize,polap,smin,smax)
SeisPatchOp(in,pind,psize,polap,smin,smax) = adjPatchOp(in,pind,psize,polap,smin,smax)

# CANDO:
Create a function  with something like fwd=1 and adj=0
PatchOp(IN, PARAM..., ::Val{1}) = fwdPatchOp(IN,PARAM...)
PatchOp(IN, PARAM..., ::Val{0}) = adjPatchOp(IN,PARAM...)

# fwd
function fwdPatchOp(IN::AbstractArray{T,N},
    psize::NTuple{N,Int},
    polap::NTuple{N,Int},
    smin::NTuple{N,Int},
    smax::NTuple{N,Int}) where {T <: Number,N}

# adj
function adjPatchOp(IN::AbstractArray{Array{T,N}},
    pind::Vector{NTuple{N,UnitRange{Int}}},
    psize::NTuple{N,Int},
    polap::NTuple{N,Int},
    smin::NTuple{N,Int},
    smax::NTuple{N,Int}) where {T <: Number, N}

###### Functions for plotting with patches ######
"""
"""
function SeisPlotPatches(args... ;
    indx::Int = 1,
    nrows::Int = 1,
    vmin="NULL",
    vmax="NULL",
    ylabel="time (samples)",
    xlabel="traces",
    fignum="patch$indx",
    xticks=[],
    yticks=[],
    xticklabel="NULL",
    yticklabel="NULL",
    ticksize=8,
    labelsize=14,
    titlesize=14,
    cmap="gray",
    pclip=90,
    fig_path=nothing)

titles = collect('a':'z');
ncols = Int.(floor(length(args)/nrows))

# Trial and error
fig_size = (3*ncols,7)

figure(fignum, figsize = fig_size)

pos=0
for (i,arg) in enumerate(args)
    pos += 1
    subplot(nrows,ncols,pos)
    to_plot = arg[indx]
    SeisPlotTX(to_plot,
        fignum=fignum,
        title="($(titles[i]))",
        titlesize=titlesize,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        xticks=xticks,
        yticks=yticks,
        xticklabels=xticklabel,
        yticklabels=yticklabel,
        ticksize=ticksize,
        pclip=pclip,
        ylabel=ylabel,
        xlabel=xlabel,
        labelsize=labelsize)
end

tight_layout();

fig_path == nothing ? nothing : savefig(fig_path)
end

"""

SeisResidualPlot(ref, args...)

Plots the difference betweein given sections with a reference section.
"""
function PatchResidualPlot(ref, args... ;
    indx::Int = 1,
    nrows::Int = 2,
    fact::Real=5,
    vmin="NULL",
    vmax="NULL",
    ylabel="time (samples)",
    xlabel="traces",
    fignum="patch$indx",
    xticks=[],
    yticks=[],
    xticklabel="NULL",
    yticklabel="NULL",
    ticksize=8,
    labelsize=14,
    titlesize=14,
    cmap="gray",
    pclip=90,
    fig_path=nothing)

titles = collect('a':'z');

ncols = length(args)

# Trial and error
fig_size = (3*ncols,7)

figure(fignum, figsize = fig_size)

pos_1=0
for (i,arg) in enumerate(args)
    pos_1 += 1
    pos_2 = pos_1 + ncols

    sect = arg[indx]
    diff = fact.*(ref[indx]-arg[indx])

    subplot(nrows,ncols,pos_1)

    SeisPlotTX(sect,
        fignum=fignum,
        title="($(titles[pos_1]))",
        titlesize=titlesize,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        xticks=xticks,
        yticks=yticks,
        xticklabels=xticklabel,
        yticklabels=yticklabel,
        ticksize=ticksize,
        pclip=pclip,
        ylabel=ylabel,
        xlabel=xlabel,
        labelsize=labelsize)
    
    subplot(nrows,ncols,pos_2)

    SeisPlotTX(diff,
        fignum=fignum,
        title="($(titles[pos_2]))",
        titlesize=titlesize,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        xticks=xticks,
        yticks=yticks,
        xticklabels=xticklabel,
        yticklabels=yticklabel,
        ticksize=ticksize,
        pclip=pclip,
        ylabel=ylabel,
        xlabel=xlabel,
        labelsize=labelsize)
end

tight_layout();

fig_path == nothing ? nothing : savefig(fig_path)
end
=#
#===========================================================================================#