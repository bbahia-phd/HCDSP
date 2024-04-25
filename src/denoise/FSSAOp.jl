import Base: vec

# Function to vectorize tuples
function vec(tuple::NTuple{N}) where N
    out = zeros(Int,length(tuple))
    for i in eachindex(out)
        out[i] = tuple[i]
    end
    return out
end

#### Complex-valued Fast SSA

## Lanczos-bsaed
"""

    OUT = fast_ssa_lanc(IN,k)

Computes the rank-k approximation of a frequency sclice via a Lanczos-based SSA.

# Arguments
-`IN::AbstractArray{T}`: Input frequency slice
-`k::Int`: Desired rank

# Expanded help

The implementation uses fast matrix-vector products based on FFTs.

"""
function fast_ssa_lanc(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = mbh_multiply(IN,x,flag="fwd");
    adj(x) = mbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, _, _ = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K))

    #Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        u = @view U[:,i]
        v = A(u,2);
        out += anti_diagonal_summation(u,v,L,K);
    end

    count = count_copy_times(vec(dims))

    return out ./ count
end

"""

    OUT = fast_ssa_qr(IN,k)

Computes the rank-k approximation of IN via a QR-based SSA.

# Arguments
-`IN::AbstractArray{T}`: Input frequency slice
-`k::Int`: Desired random projections.

# Expanded help

The implementation uses fast matrix-vector products based on FFTs.
"""
function fast_ssa_qr(IN::AbstractArray{T},k) where T

    # dimensions of array
    dims = size(IN);

    # OUT
    OUT = zero(IN)

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    # Random matrix for projection
    Ω = rand(prod(K),k);

    # projection
    P = zeros(T,(prod(L),k))
    for i = 1:k
        P[:,i] .= mbh_multiply(IN,Ω[:,i],flag="fwd")
    end

    # QR
    Qr,_ = qr(P)

    for i in 1:k

        q = Qr[:,i];

        t = mbh_multiply(IN,q,flag="adj");

        OUT .+= anti_diagonal_summation(q,t,L,K);        
    end

    count = count_copy_times(vec(dims))

    OUT .= OUT ./ count;

    return OUT

end


#### Quaternion-valued Fast SSA (QSSA)
## Lanczos-bsaed
"""

    OUT = fast_qssa_lanc(IN,k)

Computes the quaternionic rank-k approximation of IN via a Lanczos-based QSSA.

# Arguments
-`IN::AbstractArray{T}`: Input frequency slice
-`k::Int`: Desired rank

# Expanded help

The implementation uses fast matrix-vector products based on QFTs.

"""
function fast_qssa_lanc(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = qmbh_multiply(IN,x,flag="fwd");
    adj(x) = qmbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    #U, Bk, V = lanbpro(A,k,m=prod(L),n=prod(K),qflag=true)
    U, _, _ = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K),qflag=true)

    #Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    # for i in 1:k
    #     out += anti_diagonal_summation(Ub[:,i],V[:,i],L,K);
    # end

    # # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        u = @view U[:,i]
        v = A(u,2);
        out += anti_diagonal_sum(u,v,L,K);
    end

    count = count_copy_times(vec(dims))

    return out ./ count

end


afwd(IN,x) = vcat( qmbh_multiply(IN,x,flag="fwd"),
                   qmbh_multiply(invi.(IN),x,flag="fwd"),
                   qmbh_multiply(invj.(IN),x,flag="fwd"),
                   qmbh_multiply(invk.(IN),x,flag="fwd") )


function aadj(IN,x)
    l = floor(Int64,length(x)/4)
    ind1 = 1:l;     xo = @view x[ind1];
    ind2 = l+1:2l;  xi = @view x[ind2];
    ind3 = 2l+1:3l; xj = @view x[ind3]; 
    ind4 = 3l+1:4l; xk = @view x[ind4];

    out = qmbh_multiply(      IN, xo, flag="adj") .+ 
          qmbh_multiply(invi.(IN),xi, flag="adj") .+
          qmbh_multiply(invj.(IN),xj, flag="adj") .+
          qmbh_multiply(invk.(IN),xk, flag="adj")

    return out     
end

"""

    OUT = fast_aqssa_lanc(IN,k)

Computes the augmented quaternionic rank-k approximation of IN via a Lanczos-based QSSA.

# Arguments
-`IN::AbstractArray{T}`: Input frequency slice
-`k::Int`: Desired rank.

# Expanded help

The implementation uses fast matrix-vector products based on QFTs.
"""
function fast_aqssa_lanc(IN,k)  

    # dimensions of array
    dims = size(IN);

    # OUT
    OUTP = zero(IN);
    OUTI = zero(IN);
    OUTJ = zero(IN);
    OUTK = zero(IN);

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    A(x,i;kwargs...) = i == 1 ? afwd(IN,x) : aadj(IN,x)

    U, Bk, V = lanbpro(A,k,m=4prod(L),n=prod(K),qflag=true)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        OUTP += anti_diagonal_sum(Ub[1:prod(L),i],V[:,i],L,K);
        OUTI += anti_diagonal_sum(Ub[prod(L)+1:2prod(L), i],V[:,i],L,K);
        OUTJ += anti_diagonal_sum(Ub[2prod(L)+1:3prod(L),i],V[:,i],L,K);
        OUTK += anti_diagonal_sum(Ub[3prod(L)+1:4prod(L),i],V[:,i],L,K);
    end

    count = count_copy_times(vec(dims))

    OUTP .= 0.25 .* (OUTP .+ invi.(OUTI) .+ invj.(OUTJ) .+ invk.(OUTK) ) ./ count;

    return OUTP 
end

"""

    OUT = fast_qssa_qr(IN,k)

Computes the quaternionic rank-k approximation of IN via a QR-based QSSA.

# Arguments
-`IN::AbstractArray{T}`: Input frequency slice
-`k::Int`: Desired random projections.

# Expanded help

The implementation uses fast matrix-vector products based on QFTs.
"""
function fast_qssa_qr(IN::AbstractArray{Quaternion{T}},k) where T

    # dimensions of array
    dims = size(IN);

    # OUT
    OUT = zero(IN)

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    # Random matrix for projection
    Ω = Quaternion.(rand(prod(K),k),rand(prod(K),k),rand(prod(K),k),rand(prod(K),k));

    # projection
    P = zeros(eltype(IN),(prod(L),k))
    for i = 1:k
        P[:,i] .= qmbh_multiply(IN,Ω[:,i],flag="fwd")
    end

    # Quaternion QR
    Qr,_ = qqr(P)

    for i in 1:k

        q = @view Qr[:,i];

        t = qmbh_multiply(IN,q,flag="adj");

        OUT .+= anti_diagonal_sum(q,t,L,K);        
    end

    count = count_copy_times(vec(dims))

    OUT .= OUT ./ count;

    return OUT

end