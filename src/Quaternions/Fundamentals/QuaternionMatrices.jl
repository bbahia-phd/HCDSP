export Quaternion

Quaternion(w::AbstractArray{T},
           x::AbstractArray{T},
           y::AbstractArray{T},
           z::AbstractArray{T} ) where T <: Real = Quaternion.(promote.(w,x,y,z)...)

Quaternion(w::AbstractArray{T} ) where T <: Real = Quaternion.(w,zero(w),zero(w),zero(w))

Quaternion(x::AbstractArray{T},
           y::AbstractArray{T},
           z::AbstractArray{T} ) where T <: Real = Quaternion.(promote.(x,y,z)...)

zerosq(dims::Tuple,T=Float64) = Quaternion(zeros(T,dims...));

function mrandq(n::Int ;T = Float64 )
    # TODO: Improve zerosq based on Base.jl
    RQ = zerosq((n,n),T)
    for i in eachindex(RQ)
        RQ[i] = randq()
    end

    return RQ
end

