################### forward and adjoint operators ##################
"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function mbh_multiply(d::Array{Tv,N}, v::Array{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number,N}

    if N <= 5
        return hankel_multiplication(d, v; flag = flag)
    else
        error("Non-supported dimension $N")
    end

end

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function hankel_multiplication(d::Vector{Tv}, v::Vector{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number}


    # dimensions of array
    dims = size(d)

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

    if flag == "fwd"

        # check the size of v
        length(v) == K1 || error(DimensionMismatch("check the size of v"))

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1])

        # Fourier transform
        fft!(v_hat)

        # element-wise multiplication
        v_hat .= v_hat .* d_hat
        ifft!(v_hat)

        # return the last L1 element
        if Tv <: Real && Ti <: Real
            return real(v_hat[K1:dims[1]])
        else
            return v_hat[K1:dims[1]]
        end

    elseif flag == "adj"

        # check the size of v
        length(v) == L1 || error(DimensionMismatch("check the size of v"))

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1])

        # Fourier transform
        fft!(v_hat)

        # conjugate property
        d_tilde = copy(d_hat)
        for i1 = 1 : dims[1]
            j1 = i1 == 1 ? 1 : dims[1]-i1+2

            d_tilde[i1] = conj(d_hat[j1])
        end

        # element-wise multiplication
        v_hat .= v_hat .* d_tilde
        ifft!(v_hat)

        # return the last K1 element
        if Tv <: Real && Ti <: Real
            return real(v_hat[L1:dims[1]])
        else
            return v_hat[L1:dims[1]]
        end

    else
        error("non-surpported operation")
    end

end

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function hankel_multiplication(d::Matrix{Tv}, v::Array{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number}

    # dimensions of array
    dims = size(d)

    # determine the output either as array or vector
    if ndims(v) > 1
       vector_flag = false
    else
       vector_flag = true
    end

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

    if flag == "fwd"

        # check the size of v
        length(v) == K1*K2 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, K1, K2)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
        fft!(v_hat)

        # element-wise multiplication
        v_hat .= v_hat .* d_hat
        ifft!(v_hat)

        # return the last L1, L2 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[K1:dims[1],K2:dims[2]])
        else
            r = v_hat[K1:dims[1],K2:dims[2]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    elseif flag == "adj"

        # check the size of v
        length(v) == L1*L2 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, L1, L2)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2])
        fft!(v_hat)

        # conjugate property
        d_tilde = copy(d_hat)
        for i2 = 1 : dims[2]
            j2 = i2 == 1 ? 1 : dims[2]-i2+2

            for i1 = 1 : dims[1]
                j1 = i1 == 1 ? 1 : dims[1]-i1+2

                d_tilde[i1,i2] = conj(d_hat[j1,j2])
            end
        end

        # element-wise multiplication
        v_hat .= v_hat .* d_tilde
        ifft!(v_hat)

        # return the last K1,K2 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[L1:dims[1],L2:dims[2]])
        else
            r = v_hat[L1:dims[1],L2:dims[2]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    else
        error("non-surpported operation")
    end

end # end function

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function hankel_multiplication(d::Array{Tv,3}, v::Array{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number}

    # dimensions of array
    dims = size(d)

    # determine the output either as array or vector
    if ndims(v) > 1
       vector_flag = false
    else
       vector_flag = true
    end

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

    if flag == "fwd"

        # check the size of v
        length(v) == K1*K2*K3 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, K1, K2, K3)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
        fft!(v_hat)

        # element-wise multiplication
        v_hat .= v_hat .* d_hat
        ifft!(v_hat)

        # return the last L1,L2,L3 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[K1:dims[1],K2:dims[2],K3:dims[3]])
        else
            r = v_hat[K1:dims[1],K2:dims[2],K3:dims[3]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    elseif flag == "adj"

        # check the size of v
        length(v) == L1*L2*L3 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, L1, L2, L3)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3])
        fft!(v_hat)

        # conjugate property
        d_tilde = copy(d_hat)
        for i3 = 1 : dims[3]
            j3 = i3 == 1 ? 1 : dims[3]-i3+2

            for i2 = 1 : dims[2]
                j2 = i2 == 1 ? 1 : dims[2]-i2+2

                for i1 = 1 : dims[1]
                    j1 = i1 == 1 ? 1 : dims[1]-i1+2

                    d_tilde[i1,i2,i3] = conj(d_hat[j1,j2,j3])
                end
            end
        end

        # element-wise multiplication
        v_hat .= v_hat .* d_tilde
        ifft!(v_hat)

        # return the last K1,K2,K3 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[L1:dims[1],L2:dims[2],L3:dims[3]])
        else
            r = v_hat[L1:dims[1],L2:dims[2],L3:dims[3]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    else
        error("non-surpported operation")
    end

end # end function

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function hankel_multiplication(d::Array{Tv,4}, v::Array{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number}

    # dimensions of array
    dims = size(d)

    # determine the output either as array or vector
    if ndims(v) > 1
       vector_flag = false
    else
       vector_flag = true
    end

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

    if flag == "fwd"

        # check the size of v
        length(v) == K1*K2*K3*K4 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, K1, K2, K3, K4)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
        fft!(v_hat)

        # element-wise multiplication
        v_hat .= v_hat .* d_hat
        ifft!(v_hat)

        # return the last L1,L2,L3,L4 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]])
        else
            r = v_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    elseif flag == "adj"

        # check the size of v
        length(v) == L1*L2*L3*L4 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, L1, L2, L3, L4)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4])
        fft!(v_hat)

        # conjugate property
        d_tilde = copy(d_hat)
        for i4 = 1 : dims[4]
            j4 = i4 == 1 ? 1 : dims[4]-i4+2

            for i3 = 1 : dims[3]
                j3 = i3 == 1 ? 1 : dims[3]-i3+2

                for i2 = 1 : dims[2]
                    j2 = i2 == 1 ? 1 : dims[2]-i2+2

                    for i1 = 1 : dims[1]
                        j1 = i1 == 1 ? 1 : dims[1]-i1+2

                        d_tilde[i1,i2,i3,i4] = conj(d_hat[j1,j2,j3,j4])
                    end
                end
            end
        end

        # element-wise multiplication
        v_hat .= v_hat .* d_tilde
        ifft!(v_hat)

        # return the last K1,K2,K3,K4 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]])
        else
            r = v_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    else
        error("non-surpported operation")
    end

end # end function

"""
   multi-level Hankel matrix or its adjoint times a vector (vectorization of a
multi-dimensional array)
"""
function hankel_multiplication(d::Array{Tv,5}, v::Array{Ti}; flag="fwd") where {Tv<:Number, Ti<:Number}

    # dimensions of array
    dims = size(d)

    # determine the output either as array or vector
    if ndims(v) > 1
       vector_flag = false
    else
       vector_flag = true
    end

    # Fourier transform of d
    d_hat = fft(d)

    # compute L, K
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
    L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

    if flag == "fwd"

        # check the size of v
        length(v) == K1*K2*K3*K4*K5 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, K1, K2, K3, K4, K5)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
        fft!(v_hat)

        # element-wise multiplication
        v_hat .= v_hat .* d_hat
        ifft!(v_hat)

        # return the last L1,L2,L3,L4,L5 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4],K5:dims[5]])
        else
            r = v_hat[K1:dims[1],K2:dims[2],K3:dims[3],K4:dims[4],K5:dims[5]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    elseif flag == "adj"

        # check the size of v
        length(v) == L1*L2*L3*L4*L5 || error(DimensionMismatch("check the size of v"))
        v         =  reshape(v, L1, L2, L3, L4, L5)

        # reverse v and padding zeros
        v_hat = reverse_order(v; n1=dims[1], n2=dims[2], n3=dims[3], n4=dims[4], n5=dims[5])
        fft!(v_hat)

        # conjugate property
        d_tilde = copy(d_hat)
        for i5 = 1 : dims[5]
            j5 = i5 == 1 ? 1 : dims[5]-i5+2

            for i4 = 1 : dims[4]
                j4 = i4 == 1 ? 1 : dims[4]-i4+2

                for i3 = 1 : dims[3]
                    j3 = i3 == 1 ? 1 : dims[3]-i3+2

                    for i2 = 1 : dims[2]
                        j2 = i2 == 1 ? 1 : dims[2]-i2+2

                        for i1 = 1 : dims[1]
                            j1 = i1 == 1 ? 1 : dims[1]-i1+2

                            d_tilde[i1,i2,i3,i4,i5] = conj(d_hat[j1,j2,j3,j4,j5])
                        end
                    end
                end
            end
        end

        # element-wise multiplication
        v_hat .= v_hat .* d_tilde
        ifft!(v_hat)

        # return the last K1,K2,K3,K4,K5 element
        if Tv <: Real && Ti <: Real
            r = real(v_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4],L5:dims[5]])
        else
            r = v_hat[L1:dims[1],L2:dims[2],L3:dims[3],L4:dims[4],L5:dims[5]]
        end

        # make output consistent with input
        if vector_flag
            return vec(r)
        else
            return r
        end

    else
        error("non-surpported operation")
    end # forward or adjoint

end # end function



"""
   Fast anti-diagonal averaging for multi-level Hankel matrix
"""
function anti_diagonal_summation(u::Vector{Tv},
                                 v::Vector{Tv},
                                 L::Tuple{Vararg{Int}},
                                 K::Tuple{Vararg{Int}}) where {Tv<:Number}

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

    # Fourier
    fft!(u_hat)
    fft!(v_hat)
    u_hat .= u_hat .* v_hat
    ifft!(u_hat)

    # convert the output to the same type as input
    if Tv <: AbstractFloat
        return real(u_hat)
    else
        return u_hat
    end

end

"""
   slow anti-diagonal summation via building multi-level hankel matrix first.
"""
function anti_diagonal_summation_slow(u::Vector{Tv}, v::Vector{Tv},
                                      L::Union{Ti,Vector{Ti}},
                                      K::Union{Ti,Vector{Ti}}) where {Tv <: Number, Ti<:Int64}

    # order of multi-level hankel matrix
    order =  length(L)
    order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

    # outer product u * v'
    H = u * v'

    if order == 1

       # check the size of input
       L[1] == length(u) || error(DimensionMismatch("check the length of u"))
       K[1] == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1

       # allocate memory for the result
       d  = zeros(Tv,N1)

       # first layer
       for j1 = 1 : K[1]
           for i1 = 1 : L[1]
               n1 = i1 + j1 - 1

               d[n1] += H[i1,j1]
           end
       end #first


    elseif order == 2

       # check the size of input
       prod(L[1:2]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:2]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2)

       # second layer
       for j2 = 1 : K[2]
           c2 = (j2-1)*K[1]
           for i2 = 1 : L[2]
               r2 = (i2-1)*L[1]
               n2 = i2 + j2 - 1

               # first layer
               for j1 = 1 : K[1]
                   c1 = c2 + j1
                   for i1 = 1 : L[1]
                       r1 = r2 + i1
                       n1 = i1 + j1 - 1

                       d[n1,n2] += H[r1,c1]
                   end
               end #first
           end
       end #second


    elseif order == 3

       # check the size of input
       prod(L[1:3]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:3]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3)

       # anti-diagonal summation
       for j3 = 1 : K[3]
           c3 = (j3-1)*K[2]*K[1]
           for i3 = 1 : L[3]
               r3 = (i3-1)*L[2]*L[1]
               n3 = i3 + j3 - 1

               # second layer
               for j2 = 1 : K[2]
                   c2 = (j2-1)*K[1]
                   for i2 = 1 : L[2]
                       r2 = (i2-1)*L[1]
                       n2 = i2 + j2 - 1

                       # first layer
                       for j1 = 1 : K[1]
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L[1]
                               r1 = r3 + r2 + i1
                               n1 = i1 + j1 - 1

                               d[n1,n2,n3] += H[r1,c1]
                           end
                       end #first
                   end
               end #second
           end
       end #thirds


    elseif order == 4

       # check the size of input
       prod(L[1:4]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:4]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4)

       # anti-diagonal summation
       for j4 = 1 : K[4]
           c4 = (j4-1)*K[3]*K[2]*K[1]
           for i4 = 1 : L[4]
               r4 = (i4-1)*L[3]*L[2]*L[1]
               n4 = i4 + j4 - 1

               # third layer
               for j3 = 1 : K[3]
                   c3 = (j3-1)*K[2]*K[1]
                   for i3 = 1 : L[3]
                       r3 = (i3-1)*L[2]*L[1]
                       n3 = i3 + j3 - 1

                       # second layer
                       for j2 = 1 : K[2]
                           c2 = (j2-1)*K[1]
                           for i2 = 1 : L[2]
                               r2 = (i2-1)*L[1]
                               n2 = i2 + j2 - 1

                               # first layer
                               for j1 = 1 : K[1]
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L[1]
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = i1 + j1 - 1

                                       d[n1,n2,n3,n4] += H[r1,c1]
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    elseif order == 5

       # check the size of input
       prod(L[1:5]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:5]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1
       N5 = L[5]+K[5]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4, N5)

       # anti-diagonal summation
       for j5 = 1 : K[5]
           c5 = (j5-1)*K[4]*K[3]*K[2]*K[1]
           for i5 = 1 : L[5]
               r5 = (i5-1)*L[4]*L[3]*L[2]*L[1]
               n5 = i5 + j5 - 1

               for j4 = 1 : K[4]
                   c4 = (j4-1)*K[3]*K[2]*K[1]
                   for i4 = 1 : L[4]
                       r4 = (i4-1)*L[3]*L[2]*L[1]
                       n4 = i4 + j4 - 1

                       # third layer
                       for j3 = 1 : K[3]
                           c3 = (j3-1)*K[2]*K[1]
                           for i3 = 1 : L[3]
                               r3 = (i3-1)*L[2]*L[1]
                               n3 = i3 + j3 - 1

                               # second layer
                               for j2 = 1 : K[2]
                                   c2 = (j2-1)*K[1]
                                   for i2 = 1 : L[2]
                                       r2 = (i2-1)*L[1]
                                       n2 = i2 + j2 - 1

                                       # first layer
                                       for j1 = 1 : K[1]
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L[1]
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = i1 + j1 - 1

                                               d[n1,n2,n3,n4,n5] += H[r1,c1]
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up-to fiveth order")
    end

    return d
end


# """
#    Fast anti-diagonal averaging for multi-level Hankel matrix
# """
# function anti_diagonal_summation(u::Vector{Tv}, v::Vector{Tv},
#                                  L::Union{Ti,Vector{Ti}}, K::Union{Ti,Vector{Ti}}
#                                 ) where {Tv<:Number, Ti<:Int64}

#     # order of hankel matrix
#     order =  length(L)
#     order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

#     if order == 1

#        # check the size of input
#        L[1] == length(u) || error(DimensionMismatch("check the length of u"))
#        K[1] == length(v) || error(DimensionMismatch("check the length of v"))

#        # the size of array
#        N1 = L[1]+K[1]-1

#        # padding zeros
#        if Tv <: AbstractFloat
#           u_hat = zeros(Complex{Tv}, N1)
#           v_hat = zeros(Complex{Tv}, N1)
#        else
#           u_hat = zeros(Tv, N1)
#           v_hat = zeros(Tv, N1)
#        end
#        u_hat[1:L[1]] .=      u
#        v_hat[1:K[1]] .= conj(v)

#        # Fourier transform
#        fft!(u_hat)
#        fft!(v_hat)
#        u_hat .= u_hat .* v_hat
#        ifft!(u_hat)


#     elseif order == 2

#        # check the size of input
#        prod(L[1:2]) == length(u) || error(DimensionMismatch("check the length of u"))
#        prod(K[1:2]) == length(v) || error(DimensionMismatch("check the length of v"))

#        # the size of array
#        N1 = L[1]+K[1]-1
#        N2 = L[2]+K[2]-1

#        # padding zeros
#        if Tv <: AbstractFloat
#           u_hat = zeros(Complex{Tv}, N1, N2)
#           v_hat = zeros(Complex{Tv}, N1, N2)
#        else
#           u_hat = zeros(Tv, N1, N2)
#           v_hat = zeros(Tv, N1, N2)
#        end
#        u_hat[1:L[1],1:L[2]] .=      reshape(u, L[1], L[2])
#        v_hat[1:K[1],1:K[2]] .= conj(reshape(v, K[1], K[2]))

#        # Fourier transform
#        fft!(u_hat)
#        fft!(v_hat)
#        u_hat .= u_hat .* v_hat
#        ifft!(u_hat)


#     elseif order == 3

#        # check the size of input
#        prod(L[1:3]) == length(u) || error(DimensionMismatch("check the length of u"))
#        prod(K[1:3]) == length(v) || error(DimensionMismatch("check the length of v"))

#        # the size of array
#        N1 = L[1]+K[1]-1
#        N2 = L[2]+K[2]-1
#        N3 = L[3]+K[3]-1

#        # padding zeros
#        if Tv <: AbstractFloat
#           u_hat = zeros(Complex{Tv}, N1, N2, N3)
#           v_hat = zeros(Complex{Tv}, N1, N2, N3)
#        else
#           u_hat = zeros(Tv, N1, N2, N3)
#           v_hat = zeros(Tv, N1, N2, N3)
#        end
#        u_hat[1:L[1],1:L[2],1:L[3]] .=      reshape(u, L[1], L[2], L[3])
#        v_hat[1:K[1],1:K[2],1:K[3]] .= conj(reshape(v, K[1], K[2], K[3]))

#        # Fourier transform
#        fft!(u_hat)
#        fft!(v_hat)
#        u_hat .= u_hat .* v_hat
#        ifft!(u_hat)


#     elseif order == 4

#        # check the size of input
#        prod(L[1:4]) == length(u) || error(DimensionMismatch("check the length of u"))
#        prod(K[1:4]) == length(v) || error(DimensionMismatch("check the length of v"))

#        # the size of array
#        N1 = L[1]+K[1]-1
#        N2 = L[2]+K[2]-1
#        N3 = L[3]+K[3]-1
#        N4 = L[4]+K[4]-1

#        # padding zeros
#        if Tv <: AbstractFloat
#           u_hat = zeros(Complex{Tv}, N1, N2, N3, N4)
#           v_hat = zeros(Complex{Tv}, N1, N2, N3, N4)
#        else
#           u_hat = zeros(Tv, N1, N2, N3, N4)
#           v_hat = zeros(Tv, N1, N2, N3, N4)
#        end
#        u_hat[1:L[1],1:L[2],1:L[3],1:L[4]] .=      reshape(u, L[1], L[2], L[3], L[4])
#        v_hat[1:K[1],1:K[2],1:K[3],1:K[4]] .= conj(reshape(v, K[1], K[2], K[3], K[4]))

#        # Fourier transform
#        fft!(u_hat)
#        fft!(v_hat)
#        u_hat .= u_hat .* v_hat
#        ifft!(u_hat)

#     elseif order == 5

#        # check the size of input
#        prod(L[1:5]) == length(u) || error(DimensionMismatch("check the length of u"))
#        prod(K[1:5]) == length(v) || error(DimensionMismatch("check the length of v"))

#        # the size of array
#        N1 = L[1]+K[1]-1
#        N2 = L[2]+K[2]-1
#        N3 = L[3]+K[3]-1
#        N4 = L[4]+K[4]-1
#        N5 = L[5]+K[5]-1

#        # padding zeros
#        if Tv <: AbstractFloat
#           u_hat = zeros(Complex{Tv}, N1, N2, N3, N4, N5)
#           v_hat = zeros(Complex{Tv}, N1, N2, N3, N4, N5)
#        else
#           u_hat = zeros(Tv, N1, N2, N3, N4, N5)
#           v_hat = zeros(Tv, N1, N2, N3, N4, N5)
#        end
#        u_hat[1:L[1],1:L[2],1:L[3],1:L[4],1:L[5]] .=      reshape(u, L[1], L[2], L[3], L[4], L[5])
#        v_hat[1:K[1],1:K[2],1:K[3],1:K[4],1:K[5]] .= conj(reshape(v, K[1], K[2], K[3], K[4], K[5]))

#        # Fourier transform
#        fft!(u_hat)
#        fft!(v_hat)
#        u_hat .= u_hat .* v_hat
#        ifft!(u_hat)

#     else
#        error("only support up-to fiveth order")
#     end

#     # convert the output to the same type as input
#     if Tv <: AbstractFloat
#        return real(u_hat)
#     else
#        return u_hat
#     end

# end