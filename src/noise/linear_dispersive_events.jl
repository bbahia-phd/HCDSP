function linear_dispersive_events(dt,nt,h,vmin,vmax,amp,fd,fmin,fmax,f0) 
    
    # Offsets and events
    nh = length(h);
    ns = length(vmin);
    
    # Padding
    nfft = 4*nextpow(2,nt);
    
    # Allocation
    OUTF = complex.(zeros(nfft,nh))
     OUT = zeros(nfft,nh)

    # Wavelet
    wav = complex.(ricker(f0=f0, dt=dt));
     nw = length(wav);
    append!(wav,zeros(nfft-nw))
    fft!(wav,1)
    
    # Freq range
    ω_range = freq_indexes(fmin,fmax,dt,nfft[1])
    
    # Imaginary unit
    im = complex(0,1)
    
    # Loop over (positive) freqs
    @inbounds for iω in ω_range
        # Freq
        f = (iω-1)/nfft/dt;
        ω = 2*pi*f;
        for i in 1:nh, j in 1:ns
            vc = vmin[j] + (vmax[j]-vmin[j]) ./ sqrt(1 .+ (f/fd[j])^4);  # Understand that this
             p = 1 ./ vc;                                                # isn't anything physical
             t = p*abs(h[i]);                                            # but rather a theoretical way to 
                                                                         # approximate dispersive events
            @inbounds OUTF[iω,i] += exp(-im*ω*t) * wav[iω] * amp[j]
        end                                    
    end
    
    # Symmetries
    conj_symmetry!(OUTF);
    
    # ifft and truncation
    OUT = real( ifft!(OUTF,1) )[1:nt,:];

    # Return
    return OUT 
end

######################################################
function irregular_hyperbolic_events(dt,nt,h,tau,vel,amp,f0,fmin,fmax)
    
    # Padding and allocation
    nfft = (2 * nextpow(2, nt));
    
    # Offsets and events
    nh = length(h)
    ns = length(vel)
    
    # Allocation
    OUT  = zeros(nt,nh)
    OUTF = complex.(zeros(nfft,nh))
    
    # Wavelet
    wav = complex.(ricker(f0=f0, dt=dt));
     nw = length(wav);
    append!(wav,zeros(nfft-nw))
    fft!(wav,1)
    
    # Freq range
    ω_range = freq_indexes(fmin,fmax,dt,nfft)
    
    # Imaginary unit
    im = complex(0,1)
    
    # Loop over (positive) freqs
    @inbounds for iω in ω_range
        # Freq
        f = (iω-1)/nfft/dt;
        ω = 2*pi*f;
        for i in 1:nh, j in 1:ns
             t = sqrt(tau[j]^2 + (h[i] / vel[j])^2); 
            @inbounds OUTF[iω,i] += exp(-im*ω*t) * wav[iω] * amp[j]
        end                                   
    end
    
    # Symmetries
    conj_symmetry!(OUTF);
    
    # ifft and truncation
    OUT = real( ifft!(OUTF,1) )[1:nt,:];

    # Return
    return OUT 
end
