"""
    SeisLinearEvents(; <keyword arguments>)
Generate up to five dimensional data consisting of linear events.
# Arguments
- `ot=0.0`: first sample for the time axis in secs.
- `dt=0.004`: time sampling interval in secs.
- `nt=500`: number of time samples.
- `ox1=0.0`: first sample for the first spatial dimension in meters.
- `dx1=10.0`: sample interval for the first spatial dimension in meters.
- `nx1=100`: number of samples for the first spatial dimension.
- `ox2=0.0`: first sample for the second spatial dimension in meters.
- `dx2=10.0`: sample interval for the second spatial dimension in meters.
- `nx2=1`: number of samples for the second spatial dimension.
- `ox3=0.0`: second sample for the third spatial dimension in meters.
- `dx3=10.0`: sample interval for the third spatial dimension in meters.
- `nx3=1`: number of samples for the third spatial dimension.
- `ox4=0.0`: third sample for the fourth spatial dimension in meters.
- `dx4=10.0`: sample interval for the fourth spatial dimension in meters.
- `nx4=1`:number of samples for the fourth spatial dimension.
- `tau=[1.0, 1.5]`: intercept traveltimes for each event.
- `p1=[0.0001,-0.0003]`: Dip of events on the first dimension
- `p2,p3,p4=[0, 0]`: Dip of events on the following dimensions
- `amp=[1.0,-1.0]`: amplitudes for each linear event.
- `f0=20.0`: central frequency of wavelet for each linear event.
# Example
```julia
julia> using SeisPlot
julia> d = SeisLinearEvents(); SeisPlotTX(d);
```
Credits: Aaron Stanton, 2015
"""
function SeisLinearEvents(; ot=0.0, dt=0.004, nt=500, ox1=0.0, dx1=10.0,
                          nx1=100, ox2=0.0, dx2=10.0, nx2=1, ox3=0.0, dx3=10.0,
                          nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[1.0,1.5],
                          p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
                          amp=[1.0,-1.0], f0=20.0)


    w = ricker(dt=dt,f0=f0);
    nf = nextpow(2,nt);
    nw = length(w);
    t_delay = (nw-1)*dt/2;
    w = vcat(w, zeros(nf-nw));
    W = fft(w);
    x1 =ox1 .+ collect(0:1:nx1-1)*dx1;
    x2 =ox2 .+ collect(0:1:nx2-1)*dx2;
    x3 =ox3 .+ collect(0:1:nx3-1)*dx3;
    x4 =ox4 .+ collect(0:1:nx4-1)*dx4;
    n_events = length(p1);
    d = zeros(Float64, nt, nx1, nx2, nx3, nx4);
    D = zeros(Complex{Float64}, nf, nx1, nx2, nx3, nx4);
    nfh = round(Int, floor(nf/2)) + 1;
    wrs = collect(0:1:nfh-1)*2*pi/(nf*dt);     # Frequency in rad/sec

    for ie = 1:n_events
     for ix4 = 1:nx4
      for ix3 = 1:nx3
       for ix2 = 1:nx2
        for ix1 = 1:nx1
          for iw = 2:nfh-1
              phase = wrs[iw]*(tau[ie] + p1[ie]*x1[ix1] + p2[ie]*x2[ix2]+
                               + p3[ie]*x3[ix3] + p4[ie]*x4[ix4] - t_delay)
              D[iw:iw, ix1,ix2,ix3,ix4] .+= W[iw]*amp[ie]*exp(-im*phase)
              D[nf-iw+2,ix1,ix2,ix3,ix4] = conj(D[iw,ix1,ix2,ix3,ix4])
          end
        end
       end
      end
     end
    end
    d = ifft(D,1);
    d = real(d[1:nt,:,:,:,:]);
    # Add extent header
    # extent = SeisMain.Extent(convert(Int32, nt), convert(Int32, nx1),
    #                         convert(Int32, nx2), convert(Int32, nx3),
    #                         convert(Int32, nx4), convert(Float32, ot),
    #                         convert(Float32, ox1), convert(Float32, ox2),
    #                         convert(Float32, ox3), convert(Float32, ox4),
    #                         convert(Float32, dt), convert(Float32, dx1),
    #                         convert(Float32, dx2), convert(Float32, dx3),
    #                         convert(Float32,dx4), "Time", "ix1", "ix2", "ix3",
    #                         "ix4", "s", "index", "index", "index", "index", "")

    if nx4 == 1 && nx3 == 1 && nx2 == 1 && nx1 == 1
        d = reshape(d, round(Int,nt))
    elseif nx4 == 1 && nx3 == 1 && nx2 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1))
    elseif nx4 == 1 && nx3 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2))
    elseif nx4 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2), round(Int, nx3))
    else
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2), round(Int, nx3),
                    round(Int, nx4))
    end
#    return d, extent
return d
end

##############################################
"""
    SeisLinearIrregularEvents(; <keyword arguments>)
Generate up to five dimensional data consisting of irregular linear events.
# Arguments
- `ot=0.0`: first sample for the time axis in secs.
- `dt=0.004`: time sampling interval in secs.
- `nt=500`: number of time samples.
- `ox1=0.0`: first sample for the first spatial dimension in meters.
- `dx1=10.0`: sample interval for the first spatial dimension in meters.
- `nx1=100`: number of samples for the first spatial dimension.
- `ox2=0.0`: first sample for the second spatial dimension in meters.
- `dx2=10.0`: sample interval for the second spatial dimension in meters.
- `nx2=1`: number of samples for the second spatial dimension.
- `ox3=0.0`: second sample for the third spatial dimension in meters.
- `dx3=10.0`: sample interval for the third spatial dimension in meters.
- `nx3=1`: number of samples for the third spatial dimension.
- `ox4=0.0`: third sample for the fourth spatial dimension in meters.
- `dx4=10.0`: sample interval for the fourth spatial dimension in meters.
- `nx4=1`:number of samples for the fourth spatial dimension.
- `tau=[1.0, 1.5]`: intercept traveltimes for each event.
- `p1=[0.0001,-0.0003]`: Dip of events on the first dimension
- `p2,p3,p4=[0, 0]`: Dip of events on the following dimensions
- `amp=[1.0,-1.0]`: amplitudes for each linear event.
- `f0=20.0`: central frequency of wavelet for each linear event.

# Example
```julia
julia> using SeisPlot
julia> d = SeisLinearIrregularEvents(); SeisPlotTX(d);
```
"""
function SeisLinearIrregularEvents(; ot=0.0, dt=0.004, nt=500,
                        x1 = collect(0:30-1) .* 10,
                        x2 = collect(0:30-1) .* 10,
                        x3 = 0, x4 = 0,
                        p1=[0.0001,-0.0003],p2=[0.,0.],p3=[0.,0],p4=[0.,0.],
                        tau=[1.0,1.5],amp=[1.0,-1.0], f0=20.0)


    w  = ricker(dt=dt,f0=f0);
    nf = nextpow(2,nt);
    nw = length(w);
    t_delay = (nw-1)*dt/2;
    w = vcat(w, zeros(nf-nw));
    W = fft(w);

    nx1 = length(x1);
    nx2 = length(x2);
    nx3 = length(x3);
    nx4 = length(x4);

    (gx,gy) = ndgrid(x1,x2) 
    
    n_events = length(p1);
    d = zeros(Float64, nt, nx1, nx2, nx3, nx4);
    D = zeros(Complex{Float64}, nf, nx1, nx2, nx3, nx4);
    nfh = round(Int, floor(nf/2)) + 1;
    wrs = collect(0:1:nfh-1)*2*pi/(nf*dt);     # Frequency in rad/sec

    for ie = 1:n_events
        for ix4 = 1:nx4
            for ix3 = 1:nx3
                for ix2 = 1:nx2
                    for ix1 = 1:nx1
                        for iw = 2:nfh-1
                            M1 = exp.(-im .* wrs[iw] * x1[ix1] * p1[ie] )
                            M2 = exp.(-im .* wrs[iw] * x2[ix2] * p2[ie] )
                            M3 = exp.(-im .* wrs[iw] * x3[ix3] * p3[ie] )
                            M4 = exp.(-im .* wrs[iw] * x4[ix4] * p4[ie] )
                            MT = exp.(-im .* wrs[iw] * (tau[ie]-t_delay))
                
                            D[iw:iw, ix1,ix2,ix3,ix4] .+= W[iw]*amp[ie]*M1*M2*M3*M4*MT
                            D[nf-iw+2,ix1,ix2,ix3,ix4] = conj(D[iw,ix1,ix2,ix3,ix4])
                        end
                    end
                end
            end
        end
    end

    d = ifft(D,1);
    d = real(d[1:nt,:,:,:,:]);

    # Add extent header
    # extent = SeisMain.Extent(convert(Int32, nt), convert(Int32, nx1),
    #                         convert(Int32, nx2), convert(Int32, nx3),
    #                         convert(Int32, nx4), convert(Float32, ot),
    #                         convert(Float32, ox1), convert(Float32, ox2),
    #                         convert(Float32, ox3), convert(Float32, ox4),
    #                         convert(Float32, dt), convert(Float32, dx1),
    #                         convert(Float32, dx2), convert(Float32, dx3),
    #                         convert(Float32,dx4), "Time", "ix1", "ix2", "ix3",
    #                         "ix4", "s", "index", "index", "index", "index", "")

    if nx4 == 1 && nx3 == 1 && nx2 == 1 && nx1 == 1
        d = reshape(d, round(Int,nt))
    elseif nx4 == 1 && nx3 == 1 && nx2 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1))
    elseif nx4 == 1 && nx3 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2))
    elseif nx4 == 1
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2), round(Int, nx3))
    else
        d = reshape(d, round(Int, nt), round(Int, nx1),
                    round(Int, nx2), round(Int, nx3),
                    round(Int, nx4))
    end
#    return d, extent
return d
end
