mbh_multiply(d::AbstractArray{Quaternion{Tv}},v::AbstractArray{Quaternion{Ti}}; flag="fwd",μ=qi,side="left") where {Tv <: Real, Ti <: Real} = qmbh_multiply(d, v; flag=flag, μ=μ, side=side) 

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function qmbh_multiply(d, v; flag="fwd", μ=qi, side="left")

    # dimensions of array
    dims = size(d)

    # determine the output either as array or vector
    vector_flag = ndims(v) > 1 ? false : true

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    if flag == "fwd"

        # check the size of v
        length(v) == prod(K) || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, K)

        # reverse v and padding zeros
        v_hat  = PadOp(reverse(v), nin = K, npad = dims, flag = "fwd")
        
        v_hat .= HCDSP.qconv(d,v_hat,μ,side)
        
        # return the last L element
        ind = map((x,y) -> x:y, K, dims);
        r = v_hat[CartesianIndices(ind)]

    elseif flag == "adj"

        # check the size of v
        length(v) == prod(L) || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, L)

        # reverse v and padding zeros
        v_hat = PadOp(reverse(v), nin = L, npad = dims, flag = "fwd")

        v_hat .= HCDSP.qconv(conj(d),v_hat,μ,side)

        # return the last K element
        ind = map((x,y) -> x:y, L, dims);
        r = v_hat[CartesianIndices(ind)]
        
    else
        error("non-surpported operation")
    end    

    # make output consistent with input
    return vector_flag ? vec(r) : r
end

"""
   Fast anti-diagonal averaging for multi-level Hankel matrix
"""
function anti_diagonal_summation(u::Vector{Quaternion{Tv}},
                                 v::Vector{Quaternion{Tv}},
                                 L::Tuple{Vararg{Int}},
                                 K::Tuple{Vararg{Int}}) where {Tv<:Real}

    # order of hankel matrix
    order =  length(L)
    order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

    # check the size of input
    prod(L) == length(u) || error(DimensionMismatch("check the length of u"))
    prod(K) == length(v) || error(DimensionMismatch("check the length of v"))

    # output size
    N = L .+ K .- 1;

    # Pad arrays
    u_hat = PadOp(reshape(u,L),nin=L,npad=N,flag="fwd");
    v_hat = PadOp(reshape(conj(v),K),nin=K,npad=N,flag="fwd");

    return HCDSP.qconv(u_hat,v_hat,qi,"left")

end