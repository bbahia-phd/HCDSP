function ndshift(IN::AbstractArray{T}, δT::Matrix{T}, dt::Real) where {T}

    # Shift traces in IN by TT
    # Padding
    nin  = size(IN);
    npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
    INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

    # Allocation
    OUT = zero(IN);
    OUTF = zero(INF);

    # Fourier in time
    fft!(INF,1);

    # Define frequency sampling and Nyquist
    df   = 1/dt/npad[1];
    inyq = Int( floor(npad[1]/2) + 1 );
        
    # Spatial indexes
    indx = CartesianIndices( npad[2:end] );
    dω = 2π*df;

    # Loop over freqs
    @inbounds for iω in 1:inyq
        # shift
        ω = (iω-1)*dω;
        shift = exp.(1im*ω .* δT[indx])
        OUTF[iω,indx] .= INF[iω,indx] .* shift
    end

    # Symmetries and ifft
    conj_symmetry!(OUTF)

    # Truncate
    OUT .= PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );

    return OUT

end

"""
    Debiasing of coefficients α for a given data d following the model
    d = fwd(α). 

"""
function debias(fwd,adj,d,α)
    



end








#=
# Shift traces by TT
# Padding
IN = copy(pgd_fkt);
nin  = size(IN);
npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

# Allocation
OUT = zero(IN);
OUTF = zero(INF);

# Fourier in time
fft!(INF,1);

# Freq range
ω_range = freq_indexes(0, 1000, dt, npad[1]);

# Spatial indexes
indx = CartesianIndices( npad[2:end] );
dw = 2π/npad[1]/dt;

# Loop over freqs
@inbounds for iω in ω_range
    # shift
    w = (iω-1)*dw;
    shift = exp.(1im*w .* dT[indx])
    OUTF[iω,indx] .= INF[iω,indx] .* shift
end

# Symmetries and ifft
conj_symmetry!(OUTF)

# Truncate
OUT .= PadOp( real( ifft!(OUTF,1) ), nin = nin, npad = npad, flag="adj" );
=#

#=
##### Undo Shift traces by TT
IN = copy(OUTP);
nin  = size(IN);
npad = (2 * nextpow(2, nin[1]), nin[2:end]...);
INF  = complex.( PadOp(IN,nin = nin, npad = npad, flag="fwd") );

# Allocation
OUT2 = zero(IN);
OUTF2 = zero(INF);

# Fourier in time
fft!(INF,1);

# Freq range
ω_range = freq_indexes(0, 150, dt, npad[1]);

# Spatial indexes
indx = CartesianIndices( npad[2:end] );
dw = 2π/npad[1]/dt;

# Loop over freqs
@inbounds for iω in ω_range
    # shift
    w = (iω-1)*dw;
    shift = exp.(-1im*w .* dT[indx])
    OUTF2[iω,indx] .= INF[iω,indx] .* shift
end

# Symmetries and ifft
conj_symmetry!(OUTF2)

# Truncate
OUT2 .= PadOp( real( ifft!(OUTF2,1) ), nin = nin, npad = npad, flag="adj" );
=#
