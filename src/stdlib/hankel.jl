##################### building hankel matrices #####################

"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function hankel_matrix(d::Array{Tv,N}) where {Tv <: Number, N}
    if N <= 5
        return build_hankel_matrix(d)
      else
         error("only surpport up to five dimension")
      end
end

"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Vector{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

    # allocate memory for Hankel matrix
    H  = Array{Tv,2}(undef, L1, K1)

    # build hankel
    for j1 = 1 : K1
        for i1 = 1 : L1
            n1 = i1+j1-1
            H[i1,j1] = d[n1]
        end
    end

    return H
end


"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Matrix{Tv}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

    # allocate memory for Hankel matrix
    H  = Array{Tv,2}(undef, L2*L1, K2*K1)

    # build Hankel
    for j2 = 1 : K2
        c2 = (j2-1)*K1
        for i2 = 1 : L2
            r2 = (i2-1)*L1
            n2 = i2 + j2 - 1

            # first layer
            for j1 = 1 : K1
                c1 = c2 + j1
                for i1 = 1 : L1
                    r1 = r2 + i1
                    n1 = i1 + j1 - 1
                    H[r1, c1] = d[n1,n2]
                end
            end #first
        end
    end #second

    return H
end


"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Array{Tv,3}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

    # allocate memory for Hankel matrix
    H  = Array{Tv,2}(undef, L3*L2*L1, K3*K2*K1)

    # build Hankel
    for j3 = 1 : K3
        c3 = (j3-1)*K2*K1
        for i3 = 1 : L3
            r3 = (i3-1)*L2*L1
            n3 = i3 + j3 - 1

            # second layer
            for j2 = 1 : K2
                c2 = (j2-1)*K1
                for i2 = 1 : L2
                    r2 = (i2-1)*L1
                    n2 = i2 + j2 - 1

                    # first layer
                    for j1 = 1 : K1
                        c1 = c3 + c2 + j1
                        for i1 = 1 : L1
                            r1 = r3 + r2 + i1
                            n1 = i1 + j1 - 1
                            H[r1, c1] = d[n1,n2,n3]
                        end
                    end #first
                end
            end #second
        end
    end #third

    return H
end


"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Array{Tv,4}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

    # allocate memory for Hankel matrix
    H  = Array{Tv,2}(undef, L4*L3*L2*L1, K4*K3*K2*K1)

    # build Hankel
    for j4 = 1 : K4
        c4 = (j4-1)*K3*K2*K1
        for i4 = 1 : L4
            r4 = (i4-1)*L3*L2*L1
            n4 = i4 + j4 - 1

            # third layer
            for j3 = 1 : K3
                c3 = (j3-1)*K2*K1
                for i3 = 1 : L3
                    r3 = (i3-1)*L2*L1
                    n3 = i3 + j3 - 1

                    # second layer
                    for j2 = 1 : K2
                        c2 = (j2-1)*K1
                        for i2 = 1 : L2
                            r2 = (i2-1)*L1
                            n2 = i2 + j2 - 1

                            # first layer
                            for j1 = 1 : K1
                                c1 = c4 + c3 + c2 + j1
                                for i1 = 1 : L1
                                    r1 = r4 + r3 + r2 + i1
                                    n1 = i1 + j1 - 1
                                    H[r1, c1] = d[n1,n2,n3,n4]
                                end
                            end #first
                        end
                    end #second
                end
            end #third
        end
    end #fourth

    return H
end


"""
   build a multi-level block Hankel matrix from a multi-dimensional array
"""
function build_hankel_matrix(d::Array{Tv,5}) where {Tv <: Number}

    # determine the dimensions of an array
    dims = size(d)

    L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
    L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
    L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
    L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
    L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

    # allocate memory for Hankel matrix
    H  = Array{Tv,2}(undef, L5*L4*L3*L2*L1, K5*K4*K3*K2*K1)

    # build Hankel
    for j5 = 1 : K5
        c5 = (j5-1)*K4*K3*K2*K1
        for i5 = 1 : L5
            r5 = (i5-1)*L4*L3*L2*L1
            n5 = i5 + j5 - 1

            for j4 = 1 : K4
                c4 = (j4-1)*K3*K2*K1
                for i4 = 1 : L4
                    r4 = (i4-1)*L3*L2*L1
                    n4 = i4 + j4 - 1

                    # third layer
                    for j3 = 1 : K3
                        c3 = (j3-1)*K2*K1
                        for i3 = 1 : L3
                            r3 = (i3-1)*L2*L1
                            n3 = i3 + j3 - 1

                            # second layer
                            for j2 = 1 : K2
                                c2 = (j2-1)*K1
                                for i2 = 1 : L2
                                    r2 = (i2-1)*L1
                                    n2 = i2 + j2 - 1

                                    # first layer
                                    for j1 = 1 : K1
                                        c1 = c5 + c4 + c3 + c2 + j1
                                        for i1 = 1 : L1
                                            r1 = r5 + r4 + r3 + r2 + i1
                                            n1 = i1 + j1 - 1
                                            H[r1, c1] = d[n1,n2,n3,n4,n5]
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

    return H
end

#########################################################################
"""
   count the times of a multi-dimensional array elements get copyed when
building a multi-level hankel matrix
"""
function count_copy_times(dims::Union{Ti,Vector{Ti}}) where {Ti<:Int64}

    # determine the dimensions of an array
    N = length(dims)

    # level 1
    if N == 1

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1

       # assign value
       for i1 = 1 : dims[1]
           if i1 <= K1
              n1 = i1
           end
           if L1 > K1 && i1 == L1
              n1 = L1-1
           end
           if i1 > L1
              n1 = dims[1]-i1+1
           end

           count_num[i1] = n1
       end

    # level 2
    elseif N == 2

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       # second dimension
       for i2 = 1 : dims[2]
           if i2 <= K2
              n2 = i2
           end
           if L2 > K2 && i2 == L2
              n2 = L2-1
           end
           if i2 > L2
              n2 = dims[2]-i2+1
           end

           # first dimension
           for i1 = 1 : dims[1]
               if i1 <= K1
                  n1 = i1
               end
               if L1 > K1 && i1 == L1
                  n1 = L1-1
               end
               if i1 > L1
                  n1 = dims[1]-i1+1
               end

               count_num[i1,i2] = n2 * n1
           end
       end


    # level 3
    elseif N == 3

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       # third dimension
       for i3 = 1 : dims[3]
           if i3 <= K3
              n3 = i3
           end
           if L3 > K3 && i3 == L3
              n3 = L3-1
           end
           if i3 > L3
              n3 = dims[3]-i3+1
           end

           # second dimension
           for i2 = 1 : dims[2]
               if i2 <= K2
                  n2 = i2
               end
               if L2 > K2 && i2 == L2
                  n2 = L2-1
               end
               if i2 > L2
                  n2 = dims[2]-i2+1
               end

               # first dimension
               for i1 = 1 : dims[1]
                   if i1 <= K1
                      n1 = i1
                   end
                   if L1 > K1 && i1 == L1
                      n1 = L1-1
                   end
                   if i1 > L1
                      n1 = dims[1]-i1+1
                   end

                   count_num[i1,i2,i3] = n3 * n2 * n1
               end
           end
       end


    # level 4
    elseif N == 4

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       # fourth dimension
       for i4 = 1 : dims[4]
           if i4 <= K4
              n4 = i4
           end
           if L4 > K4 && i4 == L4
              n4 = L4-1
           end
           if i4 > L4
              n4 = dims[4]-i4+1
           end

           # third dimension
           for i3 = 1 : dims[3]
               if i3 <= K3
                  n3 = i3
               end
               if L3 > K3 && i3 == L3
                  n3 = L3-1
               end
               if i3 > L3
                  n3 = dims[3]-i3+1
               end

               # second dimension
               for i2 = 1 : dims[2]
                   if i2 <= K2
                      n2 = i2
                   end
                   if L2 > K2 && i2 == L2
                      n2 = L2-1
                   end
                   if i2 > L2
                      n2 = dims[2]-i2+1
                   end

                   # first dimension
                   for i1 = 1 : dims[1]
                       if i1 <= K1
                          n1 = i1
                       end
                       if L1 > K1 && i1 == L1
                          n1 = L1-1
                       end
                       if i1 > L1
                          n1 = dims[1]-i1+1
                       end

                       count_num[i1,i2,i3,i4] = n4 * n3 * n2 * n1
                   end
               end
           end
       end


    # level 5
    elseif N == 5

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4], dims[5])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       # assign value
       for i5 = 1 : dims[5]
           if i5 <= K5
              n5 = i5
           end
           if L5 > K5 && i5 == L5
              n5 = L5-1
           end
           if i5 > L5
              n5 = dims[5]-i5+1
           end

           # fourth dimension
           for i4 = 1 : dims[4]
               if i4 <= K4
                  n4 = i4
               end
               if L4 > K4 && i4 == L4
                  n4 = L4-1
               end
               if i4 > L4
                  n4 = dims[4]-i4+1
               end

               # third dimension
               for i3 = 1 : dims[3]
                   if i3 <= K3
                      n3 = i3
                   end
                   if L3 > K3 && i3 == L3
                      n3 = L3-1
                   end
                   if i3 > L3
                      n3 = dims[3]-i3+1
                   end

                   # second dimension
                   for i2 = 1 : dims[2]
                       if i2 <= K2
                          n2 = i2
                       end
                       if L2 > K2 && i2 == L2
                          n2 = L2-1
                       end
                       if i2 > L2
                          n2 = dims[2]-i2+1
                       end

                       # first dimension
                       for i1 = 1 : dims[1]
                           if i1 <= K1
                              n1 = i1
                           end
                           if L1 > K1 && i1 == L1
                              n1 = L1-1
                           end
                           if i1 > L1
                              n1 = dims[1]-i1+1
                           end

                           count_num[i1,i2,i3,i4,i5] = n5 * n4 * n3 * n2 * n1
                       end
                   end
               end
           end
       end


    else
       error("only support up to 5 dimensions")
    end

    return count_num
end

#########################################################################
"""
   reverse the order of the element of a multi-dimensional array and padding zeros
to the resultant multi-dimensional array
"""
function reverse_order(d::AbstractArray{Tv}; n1=0, n2=0, n3=0, n4=0, n5=0)  where {Tv <: Number}

    # get the dimensions of input array
    dims = size(d)
    N    = length(dims)

    # determine single or double precision
    if Tv <: Union{Int32, Float32, Complex{Float32}}
       Te =  Complex{Float32}
    else
       Te =  Complex{Float64}
    end

    if N == 1

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1

       # allocate complex array
       r = zeros(Te, n1)

       # first dimension
       for i1 = 1 : dims[1]
           j1 = dims[1]-i1+1

           # assign value
           r[i1] = d[j1]
       end


    elseif N == 2

      # padding zeros if needed
      n1 = n1 < dims[1] ? dims[1] : n1
      n2 = n2 < dims[2] ? dims[2] : n2

      # allocate complex array
      r = zeros(Te, n1, n2)

      # second dimension
      for i2 = 1 : dims[2]
          j2 = dims[2]-i2+1

          # first dimension
          for i1 = 1 : dims[1]
              j1 = dims[1]-i1+1

              # assign value
              r[i1,i2] = d[j1,j2]
          end
      end


    elseif N == 3

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3

       # allocate complex array
       r = zeros(Te, n1, n2, n3)

       # third dimension
       for i3 = 1 : dims[3]
           j3 = dims[3]-i3+1

           # second dimension
           for i2 = 1 : dims[2]
               j2 = dims[2]-i2+1

               # first dimension
               for i1 = 1 : dims[1]
                   j1 = dims[1]-i1+1

                   # assign value
                   r[i1,i2,i3] = d[j1,j2,j3]
               end
           end
       end


    elseif N == 4

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4

       # allocate complex array
       r = zeros(Te, n1, n2, n3, n4)

       # fourth dimension
       for i4 = 1 : dims[4]
           j4 = dims[4]-i4+1

           # third dimension
           for i3 = 1 : dims[3]
               j3 = dims[3]-i3+1

               # second dimension
               for i2 = 1 : dims[2]
                   j2 = dims[2]-i2+1

                   # first dimension
                   for i1 = 1 : dims[1]
                       j1 = dims[1]-i1+1

                       # assign value
                       r[i1,i2,i3,i4] = d[j1,j2,j3,j4]
                   end
               end
           end
       end


    elseif N == 5

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4
       n5 = n5 < dims[5] ? dims[5] : n5

       # allocate complex array
       r = zeros(Te, n1, n2, n3, n4, n5)

       # fiveth dimension
       for i5 = 1 : dims[5]
           j5 = dims[5]-i5+1

           # fourth dimension
           for i4 = 1 : dims[4]
               j4 = dims[4]-i4+1

               # third dimension
               for i3 = 1 : dims[3]
                   j3 = dims[3]-i3+1

                   # second dimension
                   for i2 = 1 : dims[2]
                       j2 = dims[2]-i2+1

                       # first dimension
                       for i1 = 1 : dims[1]
                           j1 = dims[1]-i1+1

                           # assign value
                           r[i1,i2,i3,i4,i5] = d[j1,j2,j3,j4,j5]
                       end
                   end
               end
           end
       end


    else
       error("only support up to 5 dimensions")
    end

    return r
end



#########################################################################
"""
   reverse the order of the element of a multi-dimensional array and padding zeros
to the resultant multi-dimensional array
"""
function conjugate_reverse_order(d::Array{Tv}; n1=0, n2=0, n3=0, n4=0, n5=0)  where {Tv <: Number}

    # get the dimensions of input array
    dims = size(d)
    N    = length(dims)

    # determine single or double precision
    if Tv <: Union{Int32, Float32, Complex{Float32}, Quaternion{Float32}}
       Te = Tv == Quaternion{Float32} ? Quaternion{Float32} : Complex{Float32}
    else
       Te = Tv == Quaternion{Float64} ? Quaternion{Float64} : Complex{Float64}
    end

    if N == 1

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1

       # allocate complex array
       r = zeros(Te, n1)

       # first dimension
       for i1 = 1 : dims[1]
        j1 = i1 == 1 ? 1 : dims[1]-i1+2

        # assign value
        r[i1] = conj(d[j1])
       end


    elseif N == 2

      # padding zeros if needed
      n1 = n1 < dims[1] ? dims[1] : n1
      n2 = n2 < dims[2] ? dims[2] : n2

      # allocate complex array
      r = zeros(Te, n1, n2)

      # second dimension
      for i2 = 1 : dims[2]
        j2 = i2 == 1 ? 1 : dims[2]-i2+2

        for i1 = 1 : dims[1]
            j1 = i1 == 1 ? 1 : dims[1]-i1+2

              # assign value
              r[i1,i2] = conj(d[j1,j2])
          end
      end


    elseif N == 3

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3

       # allocate complex array
       r = zeros(Te, n1, n2, n3)

       for i3 = 1 : dims[3]
        j3 = i3 == 1 ? 1 : dims[3]-i3+2

        for i2 = 1 : dims[2]
            j2 = i2 == 1 ? 1 : dims[2]-i2+2

            for i1 = 1 : dims[1]
                j1 = i1 == 1 ? 1 : dims[1]-i1+2

                # assign value
                r[i1,i2,i3] = conj(d[j1,j2,j3])
               end
           end
       end


    elseif N == 4

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4

       # allocate complex array
       r = zeros(Te, n1, n2, n3, n4)

       for i4 = 1 : dims[4]
        j4 = i4 == 1 ? 1 : dims[4]-i4+2

        for i3 = 1 : dims[3]
            j3 = i3 == 1 ? 1 : dims[3]-i3+2

            for i2 = 1 : dims[2]
                j2 = i2 == 1 ? 1 : dims[2]-i2+2

                for i1 = 1 : dims[1]
                    j1 = i1 == 1 ? 1 : dims[1]-i1+2

                       # assign value
                       r[i1,i2,i3,i4] = conj(d[j1,j2,j3,j4])
                   end
               end
           end
       end


    elseif N == 5

       # padding zeros if needed
       n1 = n1 < dims[1] ? dims[1] : n1
       n2 = n2 < dims[2] ? dims[2] : n2
       n3 = n3 < dims[3] ? dims[3] : n3
       n4 = n4 < dims[4] ? dims[4] : n4
       n5 = n5 < dims[5] ? dims[5] : n5

       # allocate complex array
       r = zeros(Te, n1, n2, n3, n4, n5)

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

                           # assign value
                           r[i1,i2,i3,i4,i5] = conj(d[j1,j2,j3,j4,j5])
                       end
                   end
               end
           end
       end


    else
       error("only support up to 5 dimensions")
    end

    return r
end
