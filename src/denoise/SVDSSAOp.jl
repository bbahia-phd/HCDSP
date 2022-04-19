#import TSVD: tsvd

# SVD-based rank-reduction
function rank_reduction(IN::AbstractArray{T},k::Int) where {T <: Number}
    U,_,_ = tsvd(IN,k)
    return U*U'*IN
end

# QSVD-based rank-reduction
function rank_reduction(IN::AbstractArray{Quaternion{T}},k::Int) where {T <: Number}
    Q = fwd_complex_adjoint(IN);
    U,_,_ = tsvd(Q,k)
    OUT = U*U'*Q
    return bck_complex_adjoint(OUT) 
end

function SVDSSAOp(IN::AbstractArray{T,N}, k::Int) where {T,N}
# IN is the frequency slice
# k is its expected rank
    
    # Hankelize
    H = HankelOp(IN);


    # Rank reduction
    H .= rank_reduction(H,k);
    
   # Averaging
   OUT = AveragingOp(H,size(IN))
   
   return OUT
end

    # H = vcat( H,
    #         invi.(H),
    #         invj.(H),
    #         invk.(H) );

    # OUT = AveragingOp(H[1:Int(size(H,1)/4),:],size(IN))

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