# Extract operators from base
import Base: +, -, *, /, ^,                           # For algebraic
    real, imag, conj, inv,                   # For math
    exp, sin, cos, log, sqrt,                # For math
    abs, abs2,                               # For math
    float, isreal, isinteger,isfinite,       # For math
    isnan, isinf, iszero, isone,             # For math   
    convert, promote_rule,                   # For code-typing 
    show                                     # For displaying
   
# Extract operators from other packages
import Random                                         # Random numbers
import LinearAlgebra: pinv,norm                       # Matrix operations

# Export a quaternion
export Quaternion, quaternion, quater, show
export qim, qjm, qkm, copy!

# Define a (full) quaternion q = qs + qi*i + qj*j + qk*k
struct Quaternion{T<:Real} <: Number
    qs::T
    qi::T
    qj::T
    qk::T
end

# Imaginary unities
const qim = Quaternion(0,1,0,0)
const qjm = Quaternion(0,0,1,0)
const qkm = Quaternion(0,0,0,1)

# Include functions
include("QuaternionOverwrite.jl")                    # Julia functions to quaternion
include("QuaternionAlgebra.jl")                      # Basic (quaternion) algebra
#include("QuaternionForms.jl")                       # Subtypes and forms (Pure, Cayley...)

# Type-constructors
Quaternion(w::Real,x::Real,y::Real,z::Real) = Quaternion(promote(w,x,y,z)...)
Quaternion(x::Real) = Quaternion(x,zero(x),zero(x),zero(x))
Quaternion(x::Complex) = Quaternion(real(x),imag(x),zero(real(x)),zero(real(x)))
Quaternion(x::Complex,y::Complex) = Quaternion(real(x),imag(x),real(y),imag(y))
Quaternion(x::T,y::T) where T <: Real = error("This is ambiguous unless inputs are Complex.")
Quaternion(x::T,y::T,z::T) where T <: Real = Quaternion(zero(x),x,y,z)

# Aliases
quaternion(x) = Quaternion(x)
quaternion(x,y) = Quaternion(x,y)
quaternion(x,y,z) = Quaternion(x,y,z)
quaternion(w,x,y,z) = Quaternion(w,x,y,z)
quaternion(x::Complex) = Quaternion(x)
quaternion(x::Quaternion) = x

quater(x) = Quaternion(x)
quater(x,y) = Quaternion(x,y)
quater(x,y,z) = Quaternion(x,y,z)
quater(w,x,y,z) = Quaternion(w,x,y,z)
quater(x::Complex) = Quaternion(x)
quater(x::Quaternion) = x

# Type-conversion
convert(::Type{Quaternion}, x::Real) = Quaternion(x)
convert(::Type{Quaternion}, x::Complex) = Quaternion(x)
convert(::Type{Quaternion{T}}, x::Real) where T <: Real = Quaternion(convert(T,x))
convert(::Type{Quaternion{T}}, x::Complex) where T <: Real = Quaternion(convert(Complex{T},x))
convert(::Type{Quaternion{T}}, q::Quaternion) where T <: Real  = Quaternion{T}(convert(T,q.qs), convert(T,q.qi), convert(T,q.qj), convert(T,q.qk))
convert(::Type{Quaternion{T}}, x::Quaternion{T}) where T <: Real = x
convert(::Type{T}, q::Quaternion) where T <: Real = (iszero(q.qi) && iszero(q.qj) && iszero(q.qk)) ? convert(T,q.qs) : throw(InexactError())
convert(::Type{T}, x::Quaternion{Tx}) where {T <: Real, Tx <: Real} = Quaternion(T.(x.qs),T.(x.qi),T.(x.qj),T.(x.qk))

# Promotion rules
promote_rule(::Type{Quaternion}, ::Type{T}) where T <: Real = Quaternion
promote_rule(::Type{Quaternion{T}}, ::Type{T}) where T <: Real = Quaternion{T}
promote_rule(::Type{Quaternion}, ::Type{Complex{T}}) where T <: Real = Quaternion
promote_rule(::Type{Quaternion{T}}, ::Type{Complex{T}}) where T <: Real = Quaternion{T}
promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = Quaternion{promote_type(T,S)}
promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T <: Real, S <: Real} = Quaternion{promote_type(T,S)}

# Functions
function show(io::IO, q::Quaternion)
    pm(x) = x < zero(x) ? " - $(-x)" : " + $(x)"
    print(io, q.qs, pm(q.qi), "i", pm(q.qj), "j", pm(q.qk), "k")
end

function quaternion(A::Array{S}, B::Array{T}, C::Array{U}, D::Array{V}) where {S <: Real, T <: Real, U <: Real, V <: Real}
    if !(size(A) == size(B) == size(C) == size(D))
        error("Arguments dimensions must match.");
    end
    
    Q = similar(A, typeof(quaternion(zero(S),zero(T),zero(U),zero(V))))
    for i in 1:length(Q)
        @inbounds Q[i] = quaternion(A[i],B[i],C[i],D[i])
    end
    
    return Q
end


function quaternion(A::Array{T}, B::Array{U}, C::Array{V}) where {T <: Real, U <: Real, V <: Real}
    if !(size(A) == size(B) == size(C))
        error("Arguments dimensions must match.");
    end
    
    Q = similar(A, typeof(quaternion(zero(T),zero(T),zero(U),zero(V))))
    for i in 1:length(Q)
        @inbounds Q[i] = quaternion(zero(A[i]),A[i],B[i],C[i])
    end
    
    return Q
end


function quaternion(A::Array{T}, B::Array{U}) where {T <: Complex, U <: Complex}
    if !(size(A) == size(B))
        error("Arguments dimensions must match.");
    end
    
    Q = similar(A, typeof(quaternion(zero(T),zero(T))))
    for i in 1:length(Q)
        @inbounds Q[i] = quaternion(A[i],B[i])
    end
    
    return Q
end         

function quaternion(A::Array{T}) where {T <: Real}
    Q = similar(A, typeof(quaternion(zero(T),zero(T),zero(T),zero(T))))
    for i in 1:length(Q)
        @inbounds Q[i] = quaternion(A[i])
    end
    return Q
end         

function quaternion(A::Array{T}) where {T <: Complex}
    Q = similar(A, typeof(quaternion(zero(T),zero(T))))
    for i in 1:length(Q)
        @inbounds Q[i] = quaternion(A[i])
    end
    return Q
end         

quaternion(x::AbstractArray{T}) where T <: Quaternion = x
quaternion(x::AbstractArray{T}) where T <: Complex = copy!(similar(x, Quaternion{real(eltype(x))}), x)
quaternion(x::AbstractArray{T}) where T <: Real    = copy!(similar(x, Quaternion{eltype(x)}), x)
