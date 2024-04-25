import FFTW: fft, ifft, fft!, ifft!

export μ0, qi, qj, qk, qfft, iqfft, fft, ifft
export Float32

const qi = normalize(quater(1,0,0));
const qj = normalize(quater(0,1,0));
const qk = normalize(quater(0,0,1));
const μ0 = quater(1,1,1) / sqrt(3)

## Multiple dispatch ##
fft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T},side,dims=1:ndims(IN)) where T <: Real = qfft(IN,convert(T,μ),side,dims)
fft(IN::AbstractArray{Quaternion{T}},dims=1:ndims(IN)) where T = qfft(IN,convert(T,qi),"left",dims)

ifft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T},side,dims=1:ndims(IN)) where T <: Real = iqfft(IN,convert(T,μ),side,dims)
ifft(IN::AbstractArray{Quaternion{T}},dims=1:ndims(IN)) where T = iqfft(IN,convert(T,qi),"left",dims)

function get_side(side::String)
    return side == "left" ? 1 : -1
end

function orthonormal_basis(μ::Quaternion)
    # Creates an orthonormal basis froma pure quaternion μ.
    # An optional η pure quaternion, parallel to μ, can also be given.

    # Define basis
    η = orthogonal(μ)
    ϵ = normalize(vector_product(μ,η))
    μ = normalize(μ)

    # The actual basis
    B = [μ.qi μ.qj μ.qk;
         η.qi η.qj η.qk;
         ϵ.qi ϵ.qj ϵ.qk]

    # Requires LinearAlgebra for that
    if norm(B'*B - I) > 1e-12
        println("Innacurate basis for QFT")
    end

    return B
end

function orthonormal_basis(μ::Quaternion,q::Quaternion)
    # Creates an orthonormal basis froma pure quaternion μ.
    # An optional η pure quaternion, parallel to μ, can also be given.

    # Define basis
    η = orthogonal(μ,q)
    ϵ = normalize(vector_product(μ,η))
    μ = normalize(μ)

    # The actual basis
    B = [μ.qi μ.qj μ.qk;
         η.qi η.qj η.qk;
         ϵ.qi ϵ.qj ϵ.qk]

    # Requires LinearAlgebra for that
    if norm(B'*B - I) > 1e-12
        println("Innacurate basis for QFT")
    end

    return B
end

function change_basis(IN::AbstractArray{T},B::AbstractArray{Tb}) where {T <: Number, Tb <: Real}
# Changes the basis of the elements in a quaternion array IN, to the pure basis μ.

    # Requires LinearAlgebra for that
    if norm(B'*B - I) > 1e-12
        println(norm(B'*B - I))
        println("Innacurate basis for QFT")
    end

    μ1 = quaternion(B[1,1],B[1,2],B[1,3])
    μ2 = quaternion(B[2,1],B[2,2],B[2,3])
    μ3 = quaternion(B[3,1],B[3,2],B[3,3])

    Pi = scalar_product.(IN,μ1)
    Pj = scalar_product.(IN,μ2)
    Pk = scalar_product.(IN,μ3)
    
    return quater.(scalar.(IN),Pi,Pj,Pk)
end

function qfft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T}=convert(T,μ0),side::String="left",dims=1:ndims(IN)) where T <: Real

    # Side
    S = get_side(side);
    # Axis
    μ = normalize(μ);
    # Get basis
    B = orthonormal_basis(μ)
    # Change basis
    INb = change_basis(IN,B)

    # Complex ffts
    C1 = fft(complex.(scalar.(INb),imagi.(INb)),dims)
    C2 = fft(complex.(imagj.(INb) , S .* imagk.(INb)),dims)
    
    # OUT back into original basis
    OUT  = quaternion(real(C1), imag(C1), real(C2), imag(C2))
    OUT .= change_basis(OUT,transpose(B))

    return OUT
end

function iqfft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T}=convert(T,μ0),side::String="left",dims=1:ndims(IN)) where T <: Real

    # Side
    S = get_side(side);
    # Axis
    μ = normalize(μ);
    # Get basis
    B = orthonormal_basis(μ)
    # Change basis
    INb = change_basis(IN,B)

    # Complex ffts
    C1 = ifft(complex.(scalar.(INb),imagi.(INb)),dims)
    C2 = ifft(complex.(imagj.(INb) , S .* imagk.(INb)),dims)
    
    # OUT back into original basis
    OUT  = quaternion(real(C1), imag(C1), real(C2), imag(C2))
    OUT .= change_basis(OUT,transpose(B))

    return OUT
end

"""

    simp,perp = symplectic_decomp(IN,μ)

Symplectic decomposition of a quaternion array.

"""
function symplectic_decomp(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T}) where T <: Real

    # Allocation
    INb = copy(IN);
    simp,perp = similar(IN),similar(IN);
    
    # Get basis
    B = orthonormal_basis(μ);

    # 
    ν = quaternion(B[2,1],B[2,2],B[2,3]);

    # actual change basis
    INb .= change_basis(INb,B);

    simp = real.(INb) + imagi.(INb)*μ
    perp = (imagj.(INb) + imagk.(INb)*μ)*ν

    return simp,perp
end


"""
    qconvl(f,g,μ)    

Left-sided quaternion spectral convolution.
"""
function qconvl(f,g,μ)
    
    # QFT
    Qf = fft(f,μ,"left");
    
    # symplectic_decomp
    fs,fp = symplectic_decomp(Qf,μ);

    # fwd qft g
    g_fwd = fft(g,μ,"left")

    # rev qft g
    g_rev = fft(g,-μ,"left")

    # convolution
    return ifft((fs .* g_fwd .+ fp .* g_rev),μ,"left")

end

"""
    qconvr(f,g,μ)    

Right-sided quaternion spectral convolution.
"""
function qconvr(f,g,μ)
    
    # QFT
    Qg = fft(g,μ,"right");
    
    # symplectic_decomp
    gs,gp = symplectic_decomp(Qg,μ);

    # fwd qft g
    f_fwd = fft(f,μ,"right")

    # rev qft g
    f_rev = fft(f,-μ,"right")

    # convolution
    return ifft((gs .* f_fwd .+ gp .* f_rev), μ, "right")

end

"""
    h = qconv(f,g,μ,side)

Spectral convolution of quaternion functions f and g.

"""
qconv(f,g,μ,side="left") = side == "left" ? qconvl(f,g,μ) : qconvr(f,g,μ)


function conv(f::Vector{Quaternion{T}},g::Vector{Quaternion{T}},μ=qi,side="left") where T

    
    lf = length(f)
    lg = length(g)

    lo = lf + lg - 1

    f_pad = PadOp(f,nin=size(f),npad=(lo,),flag="fwd");
    g_pad = PadOp(g,nin=size(f),npad=(lo,),flag="fwd");

    return qconv(f_pad,g_pad,μ,side)

end
