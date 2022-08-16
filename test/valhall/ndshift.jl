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


########################################################################
@everywhere function update_weights(W,r,pvals,ii; γ=2.0)
    # The function riht! multiplies W element by element.
    ε = 1e-4;
   
    p = pvals

    if p == 1.0
        @inbounds for i in eachindex(W)
            W[i] = 1 / ( abs(r[i]) +  ε);
        end

    elseif p == 2.0
        @inbounds for i in eachindex(W)
            W[i] = one(r[i]);
        end
    else
        @inbounds for i in eachindex(W)
            W[i] = 1 / (abs(r[i])^(2.0-p) +  ε)
        end
        #γ2=γ^2
        # @inbounds for i in eachindex(W)
        #     W[i] = γ2 / ( γ2 +  r[i]^2.0);
        # end        
    end
end

########################################################################
# robust thresholding
@everywhere function fk_thresh(IN::AbstractArray,sched::AbstractArray, p)

    out = copy(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad & Crop
    fwdPad(x) = PadOp(x; nin=n, npad=npad, flag="fwd");
    adjPad(x) = PadOp(x; nin=n, npad=npad, flag="adj");
    
    # Overall fwd and adj operators with transform
    FwdOp(s) = adjPad(real(ifft(s)));
    AdjOp(s) = fft(fwdPad(s))  ./ prod(npad);    

    # iht step-size
    αi = Float32(0.1);

    # iht tolerance
    εi = Float32(1e-16);    

    # Robust Iterative Hard Thresholding
    tmp,_ = riht!(FwdOp, AdjOp, out, zeros(ComplexF64,npad),
                  sched, update_weights, p;
                  α = αi,
                  maxIter=K,
                  verbose=true)

    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

########################################################################
# non robust thresholding
@everywhere function fk_thresh(IN::AbstractArray,sched::Real)

    out = similar(IN)
    n = size(IN)
    npad = 2 .* nextpow.(2,n)

    # Pad
    tmp = complex.(PadOp(IN; nin=n, npad=npad, flag="fwd"))
    
    # fft
    fft!(tmp)

    # threshold
    threshold!(tmp,sched)
    
    # Truncate
    out .= PadOp(real( ifft!(tmp) ); nin=n, npad=npad, flag="adj")
    
    # Return
    return out
end

########################################################################
# define a patching operator
function proj!(state, (psize, polap, smin, smax, sched))

    # output allocation
    out = similar(state.x);
    
    # get iteration:
    it = state.it;

    # apply patching on input
    patches,pid = fwdPatchOp(state.x, psize, polap, smin, smax);
    
    # define the fkt function with given thresholding                                                                                  
    fkt(δ) = fk_thresh(δ,sched[it])

    patches .= pmap(fkt,patches);
        
    # rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);

    # projection
    return out
end

########################################################################
# define a patching operator
function rproj!(state, (psize, polap, smin, smax, sched, pvals))

    # output allocation
    out = similar(state.x);
    
    # get iteration:
    it = state.it;

    # p value
    p = pvals[it];

    # set internal sched for RIHT
    new_sched = _schedule(sched[1], sched[it], K, "exp")
    
    # apply patching on input
    patches,pid = fwdPatchOp(state.x, psize, polap, smin, smax);
    
    # fk_thresh all patches
    if p != 2.0
        patches .= pmap(fk_thresh,
                   patches,
                   repeat([new_sched], length(patches)),
                   repeat([p],length(patches)));
    else
        # define the fkt function with given thresholding                                  
        fkt(δ) = fk_thresh(δ,sched[it])
        patches .= pmap(fkt,patches);
    end
        
    # rewrite the solution
    out .= adjPatchOp(patches, pid, psize, polap, smin, smax);
        
    # projection
    return out
end


########################################################################
# define a debias function
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
