export norm, normalize, isunit,
    involution, invi, invj, invk,
    exp, log, sin, cos

# Prepare for algebra
(+)(q::Quaternion) = Quaternion(q.qs, q.qi, q.qj, q.qk)
(-)(q::Quaternion) = Quaternion(-q.qs, -q.qi, -q.qj, -q.qk)
(*)(q::Quaternion) = error("Invalid use of operator")
(/)(q::Quaternion) = error("Invalid use of operator")

# Operations with real scalars
(+)(q::Quaternion,x::Real) = Quaternion(q.qs+x, q.qi, q.qj, q.qk)
(-)(q::Quaternion,x::Real) = Quaternion(q.qs-x, q.qi, q.qj, q.qk)
(*)(q::Quaternion,x::Real) = Quaternion(q.qs*x, q.qi*x, q.qj*x, q.qk*x)
(/)(q::Quaternion,x::Real) = Quaternion(q.qs/x, q.qi/x, q.qj/x, q.qk/x)

# Functions 
conj(q::Quaternion) = Quaternion(q.qs, -q.qi, -q.qj, -q.qk)
invi(q::Quaternion) = Quaternion(q.qs, q.qi, -q.qj, -q.qk)
invj(q::Quaternion) = Quaternion(q.qs, -q.qi, q.qj, -q.qk)
invk(q::Quaternion) = Quaternion(q.qs, -q.qi, -q.qj, q.qk)

function involution(q::Quaternion,mu::Quaternion)::Quaternion
    if !ispure(mu)
        error("Axis must be a pure quaternion")
    end
    
    mu = isunit(mu) ? mu : normalize(mu)
    
    return -mu*q*mu
end

abs2(q::Quaternion) = (q.qs*q.qs + q.qi*q.qi + q.qj*q.qj + q.qk*q.qk)
absv(q::Quaternion) = sqrt(q.qi*q.qi + q.qj*q.qj + q.qk*q.qk)    
abs(q::Quaternion) = sqrt(q.qs*q.qs + q.qi*q.qi + q.qj*q.qj + q.qk*q.qk)
inv(q::Quaternion) = conj(q)/abs2(q)

norm(q::Quaternion) = abs(q)
normalize(q::Quaternion)::Quaternion = abs(q) > 0 ? q/abs(q) : quater(0.0,1.0,0.0,0.0)
isunit(q::Quaternion)::Bool = norm(q) == 1

# Quaternion algebra
(+)(p::Quaternion,q::Quaternion) = Quaternion(p.qs+q.qs,
                                              p.qi+q.qi,
                                              p.qj+q.qj,
                                              p.qk+q.qk)

# Could be: (-)(p::Quaternion,q::Quaternion) = +(p,-q)
(-)(p::Quaternion,q::Quaternion) = Quaternion(p.qs-q.qs,
                                              p.qi-q.qi,
                                              p.qj-q.qj,
                                              p.qk-q.qk)


# Quaternion multiplication
(*)(p::Quaternion,q::Quaternion) = Quaternion(p.qs*q.qs - p.qi*q.qi - p.qj*q.qj - p.qk*q.qk,
                                              p.qs*q.qi + p.qi*q.qs + p.qj*q.qk - p.qk*q.qj,
                                              p.qs*q.qj - p.qi*q.qk + p.qj*q.qs + p.qk*q.qi,
                                              p.qs*q.qk + p.qi*q.qj - p.qj*q.qi + p.qk*q.qs)

#   (⊗)(p::Quaternion,q::Quaternion) = Quaternion(p.qs*q.qs,
#                                                  p.qi*q.qi,
#                                                  p.qj*q.qj,
#                                                  p.qk*q.qk) 

# Division (@test p*q != q*p)
(/)(p::Quaternion,q::Quaternion) =  p*inv(q)

# Exponential of a quaternion
function exp(q::Quaternion)

    exps = exp(scalar(q))
    absq = absv(q)
    fctr = absq > 0 ? exps*sin(absq)/absq : exps
    
    return Quaternion( exps*cos(absq), fctr*q.qi, fctr*q.qj, fctr*q.qk )
end

# Logarithm of a quaternion
function log(q::Quaternion)
    normq = norm(q)
    uniq = normalize(q)
    scal = scalar(uniq)
    
    mu = normalize(imag(uniq))
    theta = scal == 0 ? 0.5*pi : atan(norm(imag(uniq)),scal)
    
    return Quaternion( log(normq),  theta*mu.qi, theta*mu.qj, theta*mu.qk )
    
end                                                               

# Useful functions
(^)(p::Quaternion,q::Quaternion) = exp(p * log(q))
sqrt(q::Quaternion) = exp(0.5 * log(q))

# Sine and cosine
function sin(q::Quaternion)
    L = normalize(vector(q))
    
    return (exp(L*q) - exp(-L*q)) / (2*L)
end

function cos(q::Quaternion)
    L = normalize(vector(q))
    
    return (exp(L*q) + exp(-L*q)) / (2)
end

# Scalar and vector products
function scalar_product(μ::Quaternion,η::Quaternion)
    if μ.qs == 0 || η.qs == 0
        return μ.qi * η.qi + μ.qj * η.qj + μ.qk * η.qk
    else
        return μ.qs * η.qs + μ.qi * η.qi + μ.qj * η.qj + μ.qk * η.qk
    end
end

function vector_product(μ::Quaternion,η::Quaternion)
    # Vector (cross) product of two pure quaternions
    qi = μ.qj * η.qk - μ.qk * η.qj
    qj = μ.qk * η.qi - μ.qi * η.qk
    qk = μ.qi * η.qj - μ.qj * η.qi

    return quater(qi,qj,qk)
end

# Parallel and orthogonal quaternions
function isparallel(μ::Quaternion,η::Quaternion)
    return abs(abs(vector_product(μ,η))) < eps()
end

function orthogonal(μ::Quaternion,η::Quaternion)
    # Creates an orthogonal unit pure quaternion perpendicular to μ.
    return normalize(vector_product(μ,η))
end

function orthogonal(μ::Quaternion)
    # Creates an orthogonal unit pure quaternion perpendicular to μ.

    # Check if μ is unitary
    !isunit(μ) ? μ = normalize(μ) : nothing

    basis = BitArray(undef, (3,1))
    for (i,q) in enumerate([qim, qjm, qkm])
        isparallel(μ,q) ? basis[i] = true : nothing
    end

    if any(basis)
        if basis[1]
            η = -qkm
        elseif basis[2]
            η = -qim
        else
            η = -qjm
        end
    else
        η = quater(1,1,1)
        if isparallel(μ,η)
            η=quater(-1,-1,0)
        end
    end
    return normalize(vector_product(μ,η))
end
