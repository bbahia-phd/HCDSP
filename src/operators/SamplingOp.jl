function Sampling(IN::Vector{T}; ε::Real = 1e-10) where {T <: Number}
	return Sampling1D(IN)
end

function Sampling(IN::AbstractArray{T,N}; ε::Real = 1e-10) where {T <: Number, N}
	return Sampling2D(IN)
end

function Sampling1D(IN::Vector{T}; ε::Real = 1e-10) where {T <: Number}

    # Input size                                                              
    nin = size(IN);

    # All spatial indexes                                                         
    all_indx = CartesianIndices( nin[2:end]  );

    # Non-zero spatial indexes (sampler)                                                        
    Γ = []

    @inbounds for i in eachindex(IN)
	σ = sqrt( abs( IN[i] )  )
        σ > ε ? push!( Γ, i ) : nothing;
    end
    
    # Define the forward sampling operator
    function fwdSampler(X,IND)
        return OUT = X[[i for i in IND]]
    end
    # Forward sampling operator (handle)
    fwd(X) = fwdSampler(X,Γ)
    
    # Define adjoint sampling operator
    function adjSampler(X,IND,INPUT)
        OUT = zero(INPUT)
        j = 0;
        @inbounds for i in IND
            j += 1;
            OUT[i] = X[j]
        end
        return OUT
    end
    # Adjoint sampling operator (handle)
    adj(X) = adjSampler(X,Γ,IN)
    
    # Return function handles
    return fwd,adj
end

function Sampling2D(IN::AbstractArray{T,2}; ε::Real = 1e-10) where {T <: Number}

    # Input size                                                              
    nin = size(IN);

    # All spatial indexes                                                         
    all_indx = CartesianIndices( nin[2:end]  );

    # Non-zero spatial indexes (sampler)                                                        
    Γ = []

    @inbounds for i in all_indx
        σ = sum( abs.(IN[:,i]).^2   )
        σ > ε ? push!( Γ, i ) : nothing;
    end

    # Define the forward sampling operator
    function fwdSampler(X,IND)
        return OUT = X[:, [i for i in IND]]
    end
    fwd(X) = fwdSampler(X,Γ)
    
    # Define adjoint sampling operator
    function adjSampler(X,IND,INPUT)
        OUT = zero(INPUT)
        j = 0;
        for i in IND
            j += 1;
            OUT[:,i] .= X[:,j]
        end
        return OUT
    end
    adj(X) = adjSampler(X,Γ,IN)
    
    # Return function handles
    return fwd,adj
end


function SamplingOp(IN::AbstractArray{T,N}; ε::Real = 1e-10) where {T <: Number,N}
    # Input size
    nin = size(IN);

    # Spatial indexes
    indx = CartesianIndices( nin[2:end]  );

    # Initiate sampler
    Γ = zeros(Int,nin)

    # Safe-guard for 1D - Not sure if this is the best to do
    @inbounds for i in eachindex(IN)
        σ = sqrt( abs( IN[i] )  )
        if σ > ε
            Γ[i] = 1
        end
    end

    return Γ
end

function SamplingMtx(IN::Array{Complex{Float64}}; cutoff::Float64 = 1e-10)
    
    # Prelim
    n = prod(size(IN))

    # Initialize sampler
    OUT = sparse(I,n,n)

    # Get indexes
    IND = findall(x -> real(abs(x)) > cutoff, IN[:])
    
    # The sampler
    OUT = OUT[IND,:]

    return OUT
end

function SamplingMtx(IN::Array{Float64}; cutoff::Float64 = 1e-10)
    
    # Prelim
    n = prod(size(IN))

    # Initialize sampler
    OUT = sparse(I,n,n)

    # Get indexes
    IND = findall(x -> abs(x) > cutoff, IN[:])
    
    # The sampler
    OUT = OUT[IND,:]

    return OUT
end
