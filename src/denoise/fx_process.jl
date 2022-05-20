function fx_irregular_process(IN::AbstractArray{T},
                              dt::Real, fmin::Real, fmax::Real, grid_irr::AbstractArray, grid_reg::AbstractArray,
                              Op::Function, args...) where {T}

    # sizes
    nout = size(grid_reg)
    
    # Padding
    nin  = size(IN);
    npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
    INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

    # Allocation
     OUT = zeros(eltype(IN), nin[1] ,nout[1],nout[2]);
    OUTF = zeros(eltype(INF),npad[1],nout[1],nout[2]);
    
    # Fourier in time
    fft!(INF,1)

    # Freq range
    ω_range = freq_indexes(fmin, fmax, dt, npad[1])

    # Spatial indexes
    indxi = CartesianIndices( npad[2:end] )    
    indxo = CartesianIndices( nout[1:2] )    

    # Loop over freqs
    @inbounds for iω in ω_range
        println("Frequency $(iω) out of $(ω_range[end])")
        # Filtering
        OUTF[iω,indxo] .= Op(INF[iω,indxi], args...)
    end

    # Symmetries and ifft
    conj_symmetry!(OUTF)

    # Truncate (it is already truncated)
    OUT .= PadOp(real( ifft!(OUTF,1) ), nin = size(OUT), npad = size(OUTF), flag="adj");

    return OUT
end   

function fx_process(IN::AbstractArray{T}, dt::Real, fmin::Real, fmax::Real, Op::Function, args...) where {T}

    # Padding
    nin  = size(IN);
    npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
    INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

    # Allocation
     OUT = zero(IN);
    OUTF = zero(INF);
    
    # Fourier in time
    fft!(INF,1)

    # Freq range
    ω_range = freq_indexes(fmin, fmax, dt, npad[1])

    # Spatial indexes
    indx = CartesianIndices( npad[2:end] )    

    # Loop over freqs
    @inbounds for iω in ω_range
        # Filtering
        OUTF[iω,indx] .= Op(INF[iω,indx], args...)
    end

    # Symmetries and ifft
    conj_symmetry!(OUTF)

    # Truncate
    OUT .= PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );

    return OUT
end   

#### Frequency-domain tools ####
function freq_indexes(fmin::Real,fmax::Real,dt::Real,n::Int)
    # Define frequency sampling and Nyquist
    df   = 1/dt/n;
    inyq = Int( floor(n/2) + 1 );
    
    # Array to output
    ifmin = Int( floor(fmin / df) + 1 );
    ifmax = Int( floor(fmax / df) + 1 );

    # Safe-guards
    if ifmin < 1;    ifmin = 1;    end
    if ifmax > inyq; ifmax = inyq; end
    
    # Output
    return ifmin:ifmax
end

function conj_symmetry!(IN)
    # Dimensions
    n = size(IN)
    
    # Indexes
    inyq = Int( round(n[1]/2) );
    indx = CartesianIndices( n[2:end] )
    
    # Honor conj. symmetries
    @inbounds for iω in inyq + 2 : n[1]
        IN[iω,indx] .= conj.( IN[n[1]-iω+2,indx] )
    end
end

## Quaternion F-X Process ##
function fx_process(IN::AbstractArray{Quaternion{T}}, dt::Real, fmin::Real, fmax::Real, Op::Function, args...) where {T}

    # Padding
    nin  = size(IN);
    nfft = 2 * nextpow(2, nin[1]);
    npad = (nfft, nin[2:end]...);
    INF  = PadOp(IN,nin = nin, npad = npad, flag="fwd");

    # Allocation
     OUT = zero(IN);
    OUTF = zero(INF);
    
    # Fourier in time
    INF .= fft(INF,1)

    # Freq range
    ω_range = freq_indexes(fmin, fmax, dt, npad[1])

    # Spatial indexes
    indx = CartesianIndices( npad[2:end] )    

    # Loop over freqs
    @inbounds for iω in ω_range
        # Filtering
        OUTF[iω,indx] .= Op(INF[iω,indx], args...)
    end


    @inbounds for iω in nfft-ω_range[end]+2:nfft-ω_range[1]+1
        # Filtering
        OUTF[iω,indx] .= Op(INF[iω,indx], args...)
    end

    # Truncate
    OUT .= PadOp( ifft(OUTF,1), nin = nin, npad = npad, flag="adj" );

    return OUT
end   
