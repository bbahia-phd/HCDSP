export qqr

"""

    u,ζ = householder_vector(a,v)

Calculates a Householder vector u, with a norm of √2, and a value ζ, from the vectors a and v. v must be real.


## Extended help

The results may be used to construct a left or right Householder matrix.

"""
function householder_vector(a::AbstractVector{Ta},v::AbstractArray{Tv}) where {Ta <: Number, Tv <: Number}


    if Tv == Quaternion
        error("v cannot be a quaternion vector")
    end

    if size(a) != size(v)
        error("Input parameters must be vectors of same size")
    end

    # allocation
    u = zero(a);
    ζ = zero(Ta);
    
    α = norm(a);

    if α == 0
        u .= a .- a

        ζ = one(Ta)
        return u,ζ
    end

    # This is a scalar
    rω = transpose(a)*v
    r = abs(rω)

    if r != 0
        ζ = -rω ./ r
    else
        ζ = 1;
    end

    μ = sqrt(α .* (α + r))

    u = (a - (ζ .* v) .* α ) ./ μ

    return u,ζ
end

function qqr(A::AbstractMatrix{Quaternion{Ta}}) where Ta

    m,n = size(A);

    N = minimum([m,n]);

    Q = Quaternion.(Matrix{Ta}(I,m,m))
    R = A;

    for j = 1:N

        s0 = R[j:m,j];

        h,ζ = householder_vector(s0,Matrix{Ta}(I,m-j+1,1)[:,1])

        T = copy(R[j:m,j:n])
        R[j:m,j:n] = (1 / ζ) .* (T - h * (T'*h)')
        
        T = copy(Q[j:m,:])
        Q[j:m,:]   = (1 / ζ) .* (T - h * (T'*h)')
    end

    return Q',R
end
