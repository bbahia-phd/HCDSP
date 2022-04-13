"""
    puck(in; <keyword arguments>)

Dip estimation by Plane Wave Destruction.

See http://sepwww.stanford.edu/data/media/public/sep//prof/pvi.pdf Chapter 4.

# Arguments
-`in::Matrix{T}`: input data
-`w1=10`: length of triangle smoothing in time
-`w2=10`: length of triangle smoothing in space
-`format="angle"`

# TODO
-Create an specific function to convert from dip to angle.
"""
function puck(in::Matrix{T}; w1::Int=10, w2::Int=10) where {T}
    
    # memory allocation
    res_out = similar(in);
    dip_out = similar(in);
    coh_out = similar(in);
    
    # Discretizations
    ∂dx  = similar(in);
    ∂dt  = similar(in);

    # update ∂dx and ∂dt with proper discretization
    _discretize_x!(∂dx,in)
    _discretize_t!(∂dt,in)

    # get operators (recycle two ∂ arrays)
    ∂dtx = ∂dt.*∂dx
    ∂dt .= ∂dt.*∂dt
    ∂dx .= ∂dx.*∂dx

    # smooth operators
    triangle_smoothing!(∂dt,  nw1=w1, nw2=w2)
    triangle_smoothing!(∂dx,  nw1=w1, nw2=w2)
    triangle_smoothing!(∂dtx, nw1=w1, nw2=w2)

    # pick on a continuum
    puck!(coh_out, dip_out, ∂dx, ∂dt, ∂dtx)

    # actual pwd
    res_out .= plane_wave_destructor(in,dip_out)

    return coh_out, dip_out, res_out
end

"""
plane_wave_destroy(a,b,c,d)

Applies `a*δ_x + p*δ_t` to a section
"""
function plane_wave_destructor(d::AbstractArray{T,N}, dip::AbstractArray{T,N}; n=size(d), a::T=one(T)) where {T,N}
    
    # output allocation
    out=zero(d)
    
    # Differencing start
    s11,s12,s21,s22 = zero(d),zero(d),zero(d),zero(d);
    s11 .= -a .- dip; s12 .=  a .- dip; s21 .= -a .+ dip; s22 .=  a .+ dip;

    @inbounds for i2 in 1:n[2]-1, i1 in 1:n[1]-1
        out[i1,i2] = d[i1,i2]*s11[i1,i2] + d[i1,i2+1]*s12[i1,i2] + d[i1+1,i2]*s21[i1,i2] + d[i1+1,i2+1]*s22[i1,i2]
    end

    # boundaries
    out[end,:] .= out[end-1,:]
    out[:,end] .= out[:,end-1]

    return out
end

"""
    puck()

Estimate `coh` (coherency) and `dip` from differential operators.

Named `puck` to honor Claerbout's implementation.
"""
function puck!(coh_out, dip_out, ∂dx, ∂dt, ∂dtx; ε = eps(eltype(coh_out)))
    @inbounds for i in eachindex(∂dt)
        if abs(∂dt[i]) > ε
            coh_out[i] =  ∂dtx[i] / sqrt(∂dt[i] * ∂dx[i])
            dip_out[i] = -∂dtx[i] / ∂dt[i]
        end
    end
end

"""
2D triangle smoothing
"""
function triangle_smoothing!(in; nw1::Int=10, nw2::Int=10, n = size(in))
    @inbounds for i2 in 1:n[2]
        go = @view in[:,i2]
        in[:,i2] .= triangle_conv(go,nw1)
    end

    @inbounds for i1 in 1:n[1]
        go = @view in[i1,:]
        in[i1,:] .= triangle_conv(go,nw2)
    end
end

"""
    triangle_conv(in,nw)

Triangle smoothing.

References:
http://sepwww.stanford.edu/sep/prof/pvi/zp/paper_html/node7.html
"""
function triangle_conv(in::AbstractVector{T}, nw::Int) where T

    # input length
    n = length(in)
    
    # output allocation
    out = similar(in)
    
    # convolve with 2 rectangular boxes
    tmp = box_conv(box_conv(in,nw),nw)

    @inbounds for i in eachindex(out)
        out[i] = tmp[i+nw-1]
    end

    @inbounds for i in 1:nw-1
        out[i]     += tmp[nw-i]
        out[n-i+1] += tmp[n+nw+i-1]
    end
    
    return out
end

"""
    box_conv(in,nw)

Rectangular smoothing.

Convolves the input signal `in` with a box function of width `nw`.

References:
https://www.palomino.ch/convolution_boxcar.html
http://sepwww.stanford.edu/sep/prof/pvi/zp/paper_html/node6.html#SECTION00121000000000000000
"""
function box_conv(in::AbstractVector{T},nw::Int) where T
    
    # input length
    n = length(in)

    # conv size
    nc = n + nw -1
    
    # allocation
    out = zeros(T,nc)
    tmp = zeros(T,nw+n)
    
    # Conv begin
    tmp[1] = copy(in[1])
    @inbounds begin
        # Cummulative sum
        for i in 2:n
            tmp[i] = tmp[i-1] .+ in[i]
        end
        # boundary extension
        for i in n+1:nc
            tmp[i] = tmp[i-1]
        end

        for i = 1:nw
            out[i] = tmp[i]
        end

        for i = nw+1:nc
            out[i] = tmp[i] .- tmp[i-nw]
        end
    end

    return out ./ nw
end

"""
Specific discretization for x
"""
function _discretize_x!(∂dx,d; n = size(d))

    @inbounds for i2 in 1:n[2]-1, i1 in 1:n[1]-1
        # Differencing star 
        ∂dx[i1,i2] =  d[i1+1,i2+1] + d[i1,i2+1] - d[i1+1,i2] - d[i1,i2]
    end

    # boundaries
    ∂dx[end,:] .= ∂dx[end-1,:]
    ∂dx[:,end] .= ∂dx[:,end-1]

    return nothing
end

"""
Specific discretization for t   
"""
function _discretize_t!(∂dt, d; n = size(d))

    @inbounds for i2 in 1:n[2]-1, i1 in 1:n[1]-1
        # Differencing star
        ∂dt[i1,i2] = d[i1+1,i2] + d[i1+1,i2+1] - d[i1,i2+1] -d[i1,i2] 
    end

    # boundaries
    ∂dt[end,:] .= ∂dt[end-1,:]
    ∂dt[:,end] .= ∂dt[:,end-1]

    return nothing
end