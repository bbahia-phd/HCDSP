"""
    Ricker(; <keyword arguments>)
Creates a Ricker wavelet.
# Arguments
* `dt=0.002`: sampling interval in secs.
* `f0=20.0`: central frequency in Hz.
# Examples
```julia
julia> w = Ricker(); plot(w);
julia> w = Ricker(dt=0.004, f0=20); plot(w);
```
# Reference
Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth
ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.
"""
function ricker(; dt=0.002, f0=20.0)

    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(-nc:nc)
    b = (pi*f0*t).^2
    w = (1 .- 2 .*b).*exp.(-b)

end
