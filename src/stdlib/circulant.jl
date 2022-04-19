"""
   build a multi-level circulant matrix up to five dimension
"""
function circulant_matrix(d::AbstractArray{Tv,N}) where {Tv <: Number, N}

    if N <= 5
      return build_circulant_matrix(d)
    else
       error("only surpport up to five dimension")
    end

end

"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Vector{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    # allocate memory for Circulant matrix
    C = Array{Tv,2}(undef, dims[1], dims[1])

    for j1 = 1 : dims[1]
        for i1 = 1 : dims[1]
            n1 = i1 - j1 + 1
            n1 = n1 >= 1 ? n1 : n1+dims[1]

            C[i1,j1] = d[n1]
        end
    end

    return C
end

"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Matrix{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    # allocate memory for Circulant matrix
    C = Array{Tv,2}(undef, prod(dims[1:2]), prod(dims[1:2]))

    # build circulant
    for j2 = 1 : dims[2]
        c2 = (j2-1)*dims[1]
        for i2 = 1 : dims[2]
            r2 = (i2-1)*dims[1]
            n2 = i2 - j2 + 1
            n2 = n2 >= 1 ? n2 : n2+dims[2]

            for j1 = 1 : dims[1]
                c1 = c2+j1
                for i1 = 1 : dims[1]
                    r1 = r2+i1
                    n1 = i1 - j1 + 1
                    n1 = n1 >= 1 ? n1 : n1+dims[1]

                    C[r1,c1] = d[n1,n2]
                end
            end
        end
    end

    return C
end

"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Array{Tv,3}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    # allocate memory for Circulant matrix
    C = Array{Tv,2}(undef, prod(dims[1:3]), prod(dims[1:3]))

    # build circulant
    for j3 = 1 : dims[3]
        c3 = (j3-1)*dims[2]*dims[1]
        for i3 = 1 : dims[3]
            r3 = (i3-1)*dims[2]*dims[1]
            n3 = i3 - j3 + 1
            n3 = n3 >= 1 ? n3 : n3+dims[3]

            for j2 = 1 : dims[2]
                c2 = (j2-1)*dims[1]
                for i2 = 1 : dims[2]
                    r2 = (i2-1)*dims[1]
                    n2 = i2 - j2 + 1
                    n2 = n2 >= 1 ? n2 : n2+dims[2]

                    for j1 = 1 : dims[1]
                        c1 = c3+c2+j1
                        for i1 = 1 : dims[1]
                            r1 = r3+r2+i1
                            n1 = i1 - j1 + 1
                            n1 = n1 >= 1 ? n1 : n1+dims[1]

                            C[r1,c1] = d[n1,n2,n3]
                        end
                    end
                end
            end
        end
    end

    return C
end

"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Array{Tv,4}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    # allocate memory for Circulant matrix
    C = Array{Tv,2}(undef, prod(dims[1:4]), prod(dims[1:4]))

    # build circulant
    for j4 = 1 : dims[4]
        c4 = (j4-1)*dims[3]*dims[2]*dims[1]
        for i4 = 1 : dims[4]
            r4 = (i4-1)*dims[3]*dims[2]*dims[1]
            n4 = i4 - j4 + 1
            n4 = n4 >= 1 ? n4 : n4+dims[4]

            for j3 = 1 : dims[3]
                c3 = (j3-1)*dims[2]*dims[1]
                for i3 = 1 : dims[3]
                    r3 = (i3-1)*dims[2]*dims[1]
                    n3 = i3 - j3 + 1
                    n3 = n3 >= 1 ? n3 : n3+dims[3]

                    for j2 = 1 : dims[2]
                        c2 = (j2-1)*dims[1]
                        for i2 = 1 : dims[2]
                            r2 = (i2-1)*dims[1]
                            n2 = i2 - j2 + 1
                            n2 = n2 >= 1 ? n2 : n2+dims[2]

                            for j1 = 1 : dims[1]
                                c1 = c4+c3+c2+j1
                                for i1 = 1 : dims[1]
                                    r1 = r4+r3+r2+i1
                                    n1 = i1 - j1 + 1
                                    n1 = n1 >= 1 ? n1 : n1+dims[1]

                                    C[r1,c1] = d[n1,n2,n3,n4]
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    return C
end

"""
   build a multi-level circulant matrix up to five dimension
"""
function build_circulant_matrix(d::Array{Tv,5}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    # allocate memory for Circulant matrix
    C = Array{Tv,2}(undef, prod(dims[1:5]), prod(dims[1:5]))

    # build circulant
    for j5 = 1 : dims[5]
        c5 = (j5-1)*dims[4]*dims[3]*dims[2]*dims[1]
        for i5 = 1 : dims[5]
            r5 = (i5-1)*dims[4]*dims[3]*dims[2]*dims[1]
            n5 = i5 - j5 + 1
            n5 = n5 >= 1 ? n5 : n5+dims[5]

            for j4 = 1 : dims[4]
                c4 = (j4-1)*dims[3]*dims[2]*dims[1]
                for i4 = 1 : dims[4]
                    r4 = (i4-1)*dims[3]*dims[2]*dims[1]
                    n4 = i4 - j4 + 1
                    n4 = n4 >= 1 ? n4 : n4+dims[4]

                    for j3 = 1 : dims[3]
                        c3 = (j3-1)*dims[2]*dims[1]
                        for i3 = 1 : dims[3]
                            r3 = (i3-1)*dims[2]*dims[1]
                            n3 = i3 - j3 + 1
                            n3 = n3 >= 1 ? n3 : n3+dims[3]

                            for j2 = 1 : dims[2]
                                c2 = (j2-1)*dims[1]
                                for i2 = 1 : dims[2]
                                    r2 = (i2-1)*dims[1]
                                    n2 = i2 - j2 + 1
                                    n2 = n2 >= 1 ? n2 : n2+dims[2]

                                    for j1 = 1 : dims[1]
                                        c1 = c5+c4+c3+c2+j1
                                        for i1 = 1 : dims[1]
                                            r1 = r5+r4+r3+r2+i1
                                            n1 = i1 - j1 + 1
                                            n1 = n1 >= 1 ? n1 : n1+dims[1]

                                            C[r1,c1] = d[n1,n2,n3,n4,n5]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return C
end






















