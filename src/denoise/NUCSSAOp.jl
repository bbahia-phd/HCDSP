import LinearAlgebra: svd, Diagonal, norm
##################################################################
# Define singular value thresholding (svt)
svt_max = (X,μ) -> begin
    U,Σ,V = svd(X)
    Σt = Diagonal(max.(Σ .- μ,0))
    return U*Σt*V'
end

svt_min = (X,μ) -> begin
    U,Σ,V = svd(X)
    Σt = Diagonal(min.(Σ,μ))
    return U*Σt*V'
end
##################################################################
# Nuclear norm-based SSA function for L_2 minimization
function NUCSSAOp(IN::AbstractArray{T,M},
                  μ::Real;
                  maxIter::Int=100,
                  skip::Int=1,
		  flag::String="primal") where {T <: Number,M}
# SSA Operator solved via ADMM
# The rank-reduction operation is done via matrix shrinkage in the SVT function.
# Syntax:
#       OUT=NUCSSAOp(IN,PARAM)
# INPUTS:
#       IN - Monochromatic frequency slice
#        μ - Trade-off parameter
#        β - ADMM parameter    Should we include it?
# OUTPUTS:
#      OUT - Filtered frequency slice
#
# Breno Bahia
##

    # Define sampling forward and adjoint
    Sf,Sa = Sampling(IN);
    
    # Call primal admm for IN
    if flag == "primal"
	OUT,_ = primal_admm_hankel(Sf,
                                   Sa,
                                   IN,
                                   μ;
                                   maxIter=maxIter,
                                   skip=skip)
    elseif flag == "dual"
	OUT,_ = dual_admm_hankel(Sf,
                                 Sa,
                                 IN,
                                 μ;
                                 maxIter=maxIter,
                                 skip=skip)
        
    else
	println("flag must be either primal or dual")	
        
    end

    # return
    return OUT
end
##################################################################
# This code performs maxIter full svds
function primal_admm_hankel(fwdOp::Function,
                            adjOp::Function,
                            IN::AbstractArray{T},
                            μ::Real;
                            β::Real=μ/(2*norm(IN,2)),
                            maxIter::Int=500,
                            skip::Int=50,
                            verbose::Bool = false) where {T <: Number}
    # Hankel matrix dimensions
    N = size(IN);
    L = Int.(floor.(N ./ 2) .+ 1);
    K = Int.(N .- L .+ 1);

    # Fwd and adj Hankel operators
    AOp(x) = AveragingOp(x,N;L=L,K=K)
    HOp(x) = HankelOp(x;N=N,L=L,K=K)

    # Define cost function (for plotting)
    f_cost(X) = 0.5 * norm( fwdOp(IN .- x) )^2 + μ * sum(svdvals( HOp(x) )) 
    misfit = zeros(maxIter+1)

    # Step-sizes
    M = minimum((prod(L),prod(K)));
    σ  = 1/(M+1/β)
    τβ = 1.61*β

    # Allocate memory (TODO: find descriptive nomenclature)
    x  = copy(IN); Hx = HOp(IN)
    xl = zero(x) ; Xl = zero(Hx)
     X = zero(Hx);  Λ = zero(Hx);
        
    # Eval first iteration
#    misfit[1] = f_cost(x)
    
    # Iterate
    for k in 1:maxIter
        # Matrix shrinkage
        X .= svt_max(Λ/β - Hx, μ/β)
        
        # Get gradient
        Xl = X - Λ/β;
        xl = AOp(Xl);
        xl .+= AOp(HOp(x))
        g = xl .+ (adjOp(fwdOp(x-IN)))/β
        
        # Update solution
        x  .-= σ .* g
        Hx .= HOp(x)
        
        # Update Lagrange multipler
        Λ .-= τβ .* (Hx .+ X)
        
        # Eval cost
#        misfit[k+1] = f_cost(Hx)

        if mod(k,skip) == 0 && verbose
#            println("Iteration number ", k, " misfit ", misfit[k+1])
            println("Iteration number ", k)            
        end
    end
    
    return x,misfit
end
##################################################################
# This code performs maxIter full svds
function dual_admm_hankel(fwdOp::Function,
                          adjOp::Function,
                          IN::AbstractArray{T},
                          μ::Real;
                          β::Real=norm(IN,2)/(8*μ),
                          maxIter=500,
                          skip=50,
                          verbose::Bool=false) where {T <: Number}
    
    # Hankel matrix dimensions
    N = size(IN);
    L = Int.(floor.(N ./ 2) .+ 1);
    K = Int.(N .- L .+ 1);
    M = minimum((prod(L),prod(K)));

    # Fwd and adj Hankel operators
    AOp(x) = AveragingOp(x,N;L=L,K=K)
    HOp(x) = HankelOp(x;N=N,L=L,K=K)
    
    # Define cost function (for plotting)
    f_cost(X) = 0.5 * norm( fwdOp(IN .- x) )^2 + μ * sum(svdvals( HOp(x) )) 
    misfit = zeros(maxIter+1)

    # Parameters (ADMM and step-sizes)
    σ  = 1/M
    τβ = 1.61*β

    # Allocate memory for
    # Updating variable (Hx is its associated Hankel mtx)
    x = copy(IN); Hx = HOp(x);
    
    # Lagrange multipliers (Hλ is its associated Hankel mtx)
    λ = zero(x);  Hλ = zero(Hx);
    γ = zero(x);  

    # Eval first iteration
#    misfit[1] = f_cost(x)

    # Iterate                                                                                                                                                                            
    for k in 1:maxIter

        # γ update              
        tmp = x .+ β .* λ;
        γ  .= adjOp(fwdOp(IN .- tmp)) ./ (1 + β) 
	# Taking adjOp() here because                                                                                                                        
	# they always show up together                                                                                                                     

        # Get gradient                                                                                                                                                               
	tmp .= γ .+ λ .+ x ./ β; # Overwrite tmp                                                                                                                                          
        G = Hλ .- HOp(σ .* tmp)

        # Matrix shrinkage (SVD)                                                                                                                                                              
        Hλ .= svt_min(G, μ)
        λ  .= AOp(Hλ)

        # Actual update                                                                                                                                                                       
        x .+= τβ .* (λ + γ)
        
        # Eval cost
#        misfit[k+1] = f_cost(Hx)

        if mod(k,skip) == 0 && verbose
#            println("Iteration number ", k, " misfit ", misfit[k+1])
            println("Iteration number ", k)
        end
    end

    return x,misfit
end


######################################################
# This code performs maxIter full svds
function dual_admm_hankel_lr(Fwd,
    Adj,
    IN,
    μ,
    β;
    maxIter=500,
    skip=50)
# Steps
τβ = 1.61*β
#σ1 = 1; # Is this an assumption?
σ2 = 1/(1+β)

# Allocate memory
x  = zero(IN);  # updating variable
Hx = HankelOp(IN);

λ  = zero(x);   # Lagrange multipliers
Hλ = zero(Hx);

γ = zero(x);    # Lagrange multiplier
tmp = zero(x);  # Helper

# Initialize misfit vector
misfit = zeros(maxIter+1)
#misfit[1] = c_func(Hx,μ)

# Iterate
for k in 1:maxIter

# γ update
tmp .= Fwd(x .+ β .* λ);
γ   .= Adj((IN .- tmp) ./ (1 + β)) # Taking Adj() here because
                       # they always show up together
# Get gradient G?
tmp .= γ .+ λ .+ x ./ β; # Overwrite tmp
G = Hλ .-HankelOp(σ2 .* tmp)

# Matrix shrinkage (SVD)
Hλ .= SVT(G, μ/β)
λ  .= AveragingOp(Hλ, size(IN)) 

# Actual update
x .+= τβ .* (λ + γ)
Hx .= HankelOp(x)

#misfit[k+1] = c_func(Hx,μ)

if mod(k,skip) == 0
println("Iteration number ", k)
end
end

return x,misfit
end