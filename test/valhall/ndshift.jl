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
    npad = 1 .* nextpow.(2,n)

    # Pad & Crop
    fwdPad(x) = PadOp(x; nin=n, npad=npad, flag="fwd");
    adjPad(x) = PadOp(x; nin=n, npad=npad, flag="adj");
    
    # Overall fwd and adj operators with transform
    FwdOp(s) = adjPad(real(ifft(s)));
    AdjOp(s) = fft(fwdPad(s))  ./ prod(npad);    

    # iht step-size
    αi = Float32(0.5);

    # iht tolerance
    εi = Float32(1e-4);    

    # Robust Iterative Hard Thresholding
    tmp,_ = riht!(FwdOp, AdjOp, out, zeros(ComplexF64,npad),
                  sched, update_weights, p;
                  α = αi,
                  maxIter=K,
                  verbose=false)

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
    new_sched = _schedule(sched[it], sched[it], K, "exp")
    
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
# mute above
"""

# Arguments
-`d` the data
-`dT` the time for each trace
-`dt` time sampling 
-`δt` time taper 
"""
function mute_above(d,dT,dt; δt=0.4)

    elt = eltype(d);
    n = size(d); nt = n[1];
    out = copy(d);
    cind = CartesianIndices(dT)
    for it in 1:nt
        t = (it-1)*dt;
        for i in cind
            tmute = dT[i];
            if (t <= tmute)
                if (t > tmute - δt )
                    out[it,i] *= elt(1 - (tmute-t)/δt)
                else
                    out[it,i] *= zero(elt)
                end
            end
        end
    end

    return out
end

                

########################################################################
# mute below
"""

# Arguments
-`d` the data
-`dT` the time for each trace
-`dt` time sampling 
-`δt` time taper 
"""
function mute_below(d,dT,dt; δt=0.4)

    elt = eltype(d);
    n = size(d); nt = n[1];
    out = copy(d);
    cind = CartesianIndices(dT)
    for it in 1:nt
        t = (it-1)*dt;
        for i in cind
            tmute = dT[i];
            if (t >= tmute)
                if (t < tmute + δt )
                    out[it,i] *= elt(1 - (tmute-t)/δt)
                else
                    out[it,i] *= zero(elt)
                end
            end
        end
    end

    return out
end


########################################################################
# define a debias function
"""
    Debiasing of coefficients α for a given data d following the model
    d = fwd(α). 

"""
function debias(fwd,adj,d,α)
###################
#### Debiasing ####
n = size(OUT); npad = 1 .* nextpow.(2,n)

# Overall fwd and adj operators with transform

# Pad & Crop
fwdPad(x) = PadOp(x; nin=n, npad=npad, flag="fwd");
adjPad(x) = PadOp(x; nin=n, npad=npad, flag="adj");

fwdFT(s) = S .* adjPad( real(ifft(M .* s)) );
adjFT(s) = M .* fft(fwdPad(S .* s))  ./ prod(npad);

fwdDb(x) = fwd(fwdFT(x));
adjDb(x) = adjFT(adj(x));

# (Full-data) Fourier coefficients
α̂ = fft(fwdPad(S .* OUT));

# Set up mask for FT coefficients
M = zeros(eltype(α̂),size(α̂));
mα = median(abs.(α̂)) .* 0.5;
gtmean(x) = x > mα
M[findall(gtmean,abs.(α̂))] .= one(eltype(M)); # NB: high-freq artifacts

# thresholded initial model
x = M .* α̂;

# model blended data from new coefficients
bb = fwdDb(x);

# residual
r = b .- bb; misfit = real(dot(r,r));

# gradient
g = adjDb(r); 
gprod = real(dot(g,g));

# conj grad
p = copy(g);

# max iter for debias step
max_iter = 10;

# mist (flag,tol...)
verbose = true; εi = 1e-6;

for i in 1:max_iter
    q = fwdDb(p);

    γ = gprod / (real(dot(q,q))+1e-10);

    # model and residual update
    x .+= γ .* p
    r .-= γ .* q
    misfit = real(dot(r,r));

    # grad
    g = adjDb(r);

    gprod_new = real(dot(g,g));
    β = gprod_new / (gprod + 1e-10);
    @show stp = abs((gprod_new - gprod)/gprod)
    if stp < εi
        # last conj grad
        p .= g .+ β .* p
        break
    else
        gprod = copy(gprod_new)
    end

    # conj grad
    p .= g .+ β .* p

    verbose ? println("Iteration $i misfit = $(misfit) and gradient = $(gprod)") : nothing
end

OUTP = fwdFT(M .* x);

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
