#function add_erratic_noise(IN::Array{Float64},dt::Float64,perc::Int)
function add_erratic_noise(IN::AbstractArray{T,N},dt::Real,perc::Int) where {T,N}

    # Preliminares
        n = size(IN)  
        t = vec((0:n[1]-1)*dt);
    ntraces = prod(n[2:end])
       nerr = round(Int,ntraces*perc/100);
    
    # Allocate
    OUT = copy(reshape(IN,n[1],ntraces))
      S = zero(OUT)

    # Indexes for erratic traces
    indx = sort(sample(1:ntraces,nerr; replace=false,ordered=false))

    # Assume 5 sinusoids per trace
    @inbounds for k in 1:nerr
        s = zeros(eltype(S),length(t),1)
        for j = 1:5
            f = rand(1:30)
            f == 1 ? f = 1 + rand(1:60) : f = f;
            s += sin.(2*pi*f*t.+randn())
        end
        S[:,indx[k]] += s;
    end
                
    OUT .+= S
    return reshape(OUT,n)
end
