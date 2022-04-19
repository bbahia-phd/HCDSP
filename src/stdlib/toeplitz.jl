"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function toeplitz_matrix(d::Array{Tv,N}) where {Tv <: Number, N}

    if N <= 5
        return build_toeplitz_matrix(d)
      else
         error("only surpport up to five dimension")
      end

end

"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Vector{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

    # allocate memory for Hankel matrix
    T  = Array{Tv,2}(undef, L1, K1)

    # build hankel
    for j1 = 1 : K1
        for i1 = 1 : L1
            n1 = K1 + i1 - j1

            T[i1,j1] = d[n1]
        end
    end

    return T
end

"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Matrix{Tv}) where {Tv <: Number}


    # determine the dimensions of an array
    dims = size(d)
     
    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

    # allocate memory for Hankel matrix
    T  = Array{Tv,2}(undef, L2*L1, K2*K1)

    # build Hankel
    for j2 = 1 : K2
        c2 = (j2-1)*K1
        for i2 = 1 : L2
            r2 = (i2-1)*L1
            n2 = K2 + i2 - j2

            # first layer
            for j1 = 1 : K1
                c1 = c2 + j1
                for i1 = 1 : L1
                    r1 = r2 + i1
                    n1 = K1 + i1 - j1

                    T[r1, c1] = d[n1,n2]
                end
            end #first
        end
    end #second

    return T

end

"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Array{Tv,3}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

    # allocate memory for Hankel matrix
    T  = Array{Tv,2}(undef, L3*L2*L1, K3*K2*K1)

    # build Hankel
    for j3 = 1 : K3
        c3 = (j3-1)*K2*K1
        for i3 = 1 : L3
            r3 = (i3-1)*L2*L1
            n3 = K3 + i3 - j3

            # second layer
            for j2 = 1 : K2
                c2 = (j2-1)*K1
                for i2 = 1 : L2
                    r2 = (i2-1)*L1
                    n2 = K2 + i2 - j2

                    # first layer
                    for j1 = 1 : K1
                        c1 = c3 + c2 + j1
                        for i1 = 1 : L1
                            r1 = r3 + r2 + i1
                            n1 = K1 + i1 - j1

                            T[r1, c1] = d[n1,n2,n3]
                        end
                    end #first
                end
            end #second
        end
    end #third

    return T

end

"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Array{Tv,4}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

    # allocate memory for Hankel matrix
    T  = Array{Tv,2}(undef, L4*L3*L2*L1, K4*K3*K2*K1)

    # build Hankel
    for j4 = 1 : K4
        c4 = (j4-1)*K3*K2*K1
        for i4 = 1 : L4
            r4 = (i4-1)*L3*L2*L1
            n4 = K4 + i4 - j4

            # third layer
            for j3 = 1 : K3
                c3 = (j3-1)*K2*K1
                for i3 = 1 : L3
                    r3 = (i3-1)*L2*L1
                    n3 = K3 + i3 - j3

                    # second layer
                    for j2 = 1 : K2
                        c2 = (j2-1)*K1
                        for i2 = 1 : L2
                            r2 = (i2-1)*L1
                            n2 = K2 + i2 - j2

                            # first layer
                            for j1 = 1 : K1
                                c1 = c4 + c3 + c2 + j1
                                for i1 = 1 : L1
                                    r1 = r4 + r3 + r2 + i1
                                    n1 = K1 + i1 - j1

                                    T[r1, c1] = d[n1,n2,n3,n4]
                                end
                            end #first
                        end
                    end #second
                end
            end #third
        end
    end #fourth

    return T

end

"""
   build a multi-level Toeplitz matrix from a multi-dimensional array
"""
function build_toeplitz_matrix(d::Array{Tv,5}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
    L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

    # allocate memory for Hankel matrix
    T  = Array{Tv,2}(undef, L5*L4*L3*L2*L1, K5*K4*K3*K2*K1)

    # build Hankel
    for j5 = 1 : K5
        c5 = (j5-1)*K4*K3*K2*K1
        for i5 = 1 : L5
            r5 = (i5-1)*L4*L3*L2*L1
            n5 = K5 + i5 - j5

            for j4 = 1 : K4
                c4 = (j4-1)*K3*K2*K1
                for i4 = 1 : L4
                    r4 = (i4-1)*L3*L2*L1
                    n4 = K4 + i4 - j4

                    # third layer
                    for j3 = 1 : K3
                        c3 = (j3-1)*K2*K1
                        for i3 = 1 : L3
                            r3 = (i3-1)*L2*L1
                            n3 = K3 + i3 - j3

                            # second layer
                            for j2 = 1 : K2
                                c2 = (j2-1)*K1
                                for i2 = 1 : L2
                                    r2 = (i2-1)*L1
                                    n2 = K2 + i2 - j2

                                    # first layer
                                    for j1 = 1 : K1
                                        c1 = c5 + c4 + c3 + c2 + j1
                                        for i1 = 1 : L1
                                            r1 = r5 + r4 + r3 + r2 + i1
                                            n1 = K1 + i1 - j1

                                            T[r1, c1] = d[n1,n2,n3,n4,n5]
                                        end
                                    end #first
                                end
                            end #second
                        end
                    end #third
                end
            end #fourth
        end
    end # fifth

    return T
end