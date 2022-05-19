function decimate_traces(IN::AbstractArray{T},perc::Real) where {T}

    # Preliminares
        n = size(IN)  
    ntraces = prod(n[2:end])
       ndec = round(Int,ntraces*perc/100);
    
    # Allocate
    OUT = copy(reshape(IN,n[1],ntraces))
    
    # Indexes for decimated traces
    indx = sort(sample(1:ntraces,ndec; replace=false,ordered=false))
    
    # Decimate
    OUT[:,indx] .= 0
    
    return reshape(OUT,n)
end

function decimate_indexes(IN::AbstractArray{T},perc::Real) where {T}

    # Preliminares
    n = size(IN)  
    ntraces = prod(n[2:end])
       ndec = round(Int,ntraces*perc/100);

    # Indexes for decimated traces
    return sort(sample(1:ntraces,ndec; replace=false,ordered=false))

end