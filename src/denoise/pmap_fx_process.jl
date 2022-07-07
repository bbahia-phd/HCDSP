function pmap_fx_process(IN::AbstractArray{T,N}, dt::Real, fmin::Real, fmax::Real, Op::Function, args...) where {T,N}

    # Padding
    nin  = size(IN);
    npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
    INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

    # Allocation
     OUT = zero(IN);
    OUTF = zero(INF);
    
    # Fourier in time
    fft!(INF,1)

    # Spatial indexes
    indx = CartesianIndices( npad[2:end] )    

    # Freq range
    # Define frequency sampling and Nyquist
    df   = 1/dt/npad[1];
    inyq = Int( floor(npad[1]/2) + 1 );
    
    # Array to output
    ifmin = Int( floor(fmin / df) + 1 );
    ifmax = Int( floor(fmax / df) + 1 );

    # Safe-guards
    if ifmin < 1;    ifmin = 1;    end
    if ifmax > inyq; ifmax = inyq; end

    ω_range =  ifmin:ifmax;
    
    nω = length(ω_range);
    outf = Vector{Array{eltype(INF),N-1}}(undef,nω);

    @inbounds for (i,iω) in enumerate(ω_range)
        outf[i] = INF[iω,indx]
    end

    # map frequency slices
    outf .= pmap(Op,outf,repeat([args...], nω));

    # change to fix symmetry here
    @inbounds for (i,iω) in enumerate(ω_range)
        # Filtering
        OUTF[iω,indx] .= outf[i]
    end

    # Symmetries and ifft
    pmap_conj_symmetry!(OUTF)

    # Truncate
    OUT .= PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );

    return OUT
end   

function pmap_conj_symmetry!(IN)
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
function pmap_fx_process(IN::AbstractArray{Quaternion{T},N}, dt::Real, fmin::Real, fmax::Real, Op::Function, args...;) where {T,N}

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
    # Define frequency sampling and Nyquist
    df   = 1/dt/npad[1];
    inyq = Int( floor(npad[1]/2) + 1 );
    
    # Array to output
    ifmin = Int( floor(fmin / df) + 1 );
    ifmax = Int( floor(fmax / df) + 1 );

    # Safe-guards
    if ifmin < 1;    ifmin = 1;    end
    if ifmax > inyq; ifmax = inyq; end

    # pos & neg frequency range
    pos_ω_range = collect(ifmin:ifmax);
    neg_ω_range = collect(nfft-ifmax+2:nfft-ifmin+1);
    ω_range = vcat(pos_ω_range,neg_ω_range)
    nω = length(ω_range);

    # Spatial indexes
    indx = CartesianIndices( npad[2:end] )
    
    # vector of frequency slices
    outf = Vector{Array{eltype(INF),N-1}}(undef,nω);
    @inbounds for (i,iω) in enumerate(ω_range)
        outf[i] = INF[iω,indx]
    end

    # map frequency slices
    outf .= pmap(Op,outf,repeat([args...], nω));

    # Loop over freqs
    @inbounds for (i,iω) in enumerate(ω_range)
        # Filtering
        OUTF[iω,indx] .= outf[i]
    end

    # Truncate
    OUT .= PadOp( ifft(OUTF,1), nin = nin, npad = npad, flag="adj" );

    return OUT
end   
