using SpecialFunctions: besseli

"""
interp_ks3d(d,grid_irr,gred_reg,N,a,flag)

# Arguments
d::AbstractMatrix{T} of size nt × ntraces containing the observed data
grid_irr::AbstractMatrix{T} of size ntraces × 2 containing the position of the observed traces
grid_reg::AbstractArray{T} of size nx × ny × 2 containing the position of the desired traces
N::Int for the Kaiser interpolator
a::Int also for the Kaiser interpolator
flag::String "fwd" maps the grids from regular to irregular or "adj" from irregular to regular
"""
function interp_ks3d(d::AbstractArray{T} ,grid_irr::AbstractArray{Tg} ,grid_reg::AbstractArray{Tg} ,N::Int,a::Int,flag::String) where {T,Tg}

    if flag == "fwd"
        return interp_ks3d_fwd(d, grid_irr, grid_reg, N, a)
    else
        return interp_ks3d_adj(d, grid_irr, grid_reg, N, a)
    end

end

"""
    "fwd" maps the grids from regular to irregular or "adj" from irregular to regular
"""
function interp_ks3d_fwd(d::AbstractArray{T} ,grid_irr::AbstractMatrix{Tg} ,grid_reg::AbstractArray{Tg} ,N::Int,a::Int) where {T,Tg}

# Get spacing from ther regular grid along x and y
dx = grid_reg[2,1,1] - grid_reg[1,1,1]
dy = grid_reg[1,2,2] - grid_reg[1,1,2]

# number of regular grid points along x and y
nx,ny,_ = size(grid_reg)

# number of traces observed in d
nk = size(grid_irr,1)

# allocate drec
drec = zeros(eltype(d),nk)

# loop over every observed position
for k in 1:nk

    # these are indexes
    ia = floor(Int,(grid_irr[k,1] - grid_reg[1,1,1]) /dx ) + 1;
    ja = floor(Int,(grid_irr[k,2] - grid_reg[1,1,2]) /dy ) + 1;

    # boundaries for ia (indexes along x) and ja (indexes along y)
    min_ia = max(ia-N,1);
    max_ia = min(ia+N,nx);

    min_ja = max(ja-N,1);
    max_ja = min(ja+N,ny);

    for ix in min_ia:max_ia
        for iy in min_ja:max_ja

            # Why t and u? 
            t = (grid_irr[k,1] - grid_reg[ix,iy,1]) / dx;
            u = (grid_irr[k,2] - grid_reg[ix,iy,2]) / dy;

            # Get weights for kaiser interpolator
            Wt = kaiser_sinc(t,N+1,a) # this is a scalar (real or compelx?)
            Wu = kaiser_sinc(u,N+1,a) # this is a scalar (real or complex?)

            # fill drec
            drec[k] += Wt*Wu*d[ix,iy]
        end
    end
end
    
return drec
end

"""
    "fwd" maps the grids from regular to irregular or "adj" from irregular to regular
"""
function interp_ks3d_adj(d::AbstractArray{T} ,grid_irr::AbstractMatrix{Tg} ,grid_reg::AbstractArray{Tg},N::Int,a::Int) where {T,Tg}

# Get spacing from ther regular grid along x and y
dx = grid_reg[2,1,1] - grid_reg[1,1,1]
dy = grid_reg[1,2,2] - grid_reg[1,1,2]

# number of regular grid points along x and y
nx,ny,_ = size(grid_reg)

# number of traces observed in d
#nt,nk = size(d)
nk = size(grid_irr,1)

# allocate drec
drec = zeros(eltype(d),nx,ny)

# loop over every observed position
for k in 1:nk

    # these are indexes
    ia = floor(Int,(grid_irr[k,1] - grid_reg[1,1,1]) /dx ) + 1;
    ja = floor(Int,(grid_irr[k,2] - grid_reg[1,1,2]) /dy ) + 1;

    # boundaries for ia (indexes along x) and ja (indexes along y)
    min_ia = max(ia-N,1);
    max_ia = min(ia+N,nx);

    min_ja = max(ja-N,1);
    max_ja = min(ja+N,ny);

    for ix in min_ia:max_ia
        for iy in min_ja:max_ja

            # Why t and u? 
            t = (grid_irr[k,1] - grid_reg[ix,iy,1]) / dx;
            u = (grid_irr[k,2] - grid_reg[ix,iy,2]) / dy;

            # Get weights for kaiser interpolator
            Wt = kaiser_sinc(t,N+1,a) # this is a scalar (real or compelx?)
            Wu = kaiser_sinc(u,N+1,a) # this is a scalar (real or complex?) real I think

            # fill drec
            drec[ix,iy] += Wt*Wu*d[k]
        end
    end
end
    
return drec

end

function kaiser_sinc(t,N,a)
    return besseli(0,a*sqrt(1-(t/N).^2))/besseli(0,a) * sinc(t)
end