export μ0, qfft, iqfft

const μ0 = quater(1,1,1) / sqrt(3)

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

function qfft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T}=μ0,side::String="left",dims=1:ndims(IN)) where T <: Real

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

function iqfft(IN::AbstractArray{Quaternion{T}},μ::Quaternion{T}=μ0,side::String="left",dims=1:ndims(IN)) where T <: Real

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
