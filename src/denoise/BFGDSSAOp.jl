##################################################################
import Statistics: median

##################################################################
# BFGD-based Robust SSA
function BFGDRSSAOp(IN::AbstractArray{T,M},
                    k::Int,
                    maxIter::Int,
                    λ::Real,
                    f_cost::String,
                    args...;
                    init::Int = 2) where {T <: Number,M}
# SSA Operator via Bifactored Gradient Descent
# The rank-reduction operation is done via bifactored regression on UV'.
# Syntax:
#       OUT=BFGDSSAOp(IN,PARAM)
# INPUTS:
#       IN - Monochromatic frequency slice
#        k - Desited rank
# OUTPUTS:
#      OUT - Filtered frequency slice
#
# Breno Bahia
##

    # Define sampling
    FwdS,AdjS = Sampling(IN);
    # TODO: Sf,Sa = Sampling(IN)
    # where IN could be Vector, Matrix, Tensor
    # (AbstractArray) but there are two different sampling functions
    # one that removes missing samples and its "put back" adjoint
    # and another that is point-wise multiplication by 1 or 0
    # which should be orthogonal.
    
    # Sample the observed data (remove zeroed sample)
    IN0 = FwdS(IN);

    # Hankel matrix dimensions
    N = size(IN);
    L = Int.(floor.(N ./ 2) .+ 1);
    K = Int.(N .- L .+ 1);
         
    # Define forward and adjoint operators
    AOp(x) = FwdS( AveragingOp(x,N,L=L,K=K) )# Forward operator S(A(X))
    HOp(x) = HankelOp( AdjS(x),N=N,L=L,K=K ) # Adjoint operator A'(S'(X))
    
    # Weighting function
    update_weights = weighting_function(f_cost)

    # Initialize initial guess
    U0,V0,η = bfgd_init(AOp,HOp,IN0,prod(L),prod(K),k; init=init)
    
    # Call bfgd for IN0
    Hk,_ = w_bfgd!(AOp,
                   HOp,
                   IN0,
                   U0,
                   V0,
                   η,
                   update_weights,
                   args...;
                   λ = λ, maxIter = maxIter,
                   verbose=false, tol=1e-5)
    
    # frequency slice (reconstruct everywhere)
    OUT = AveragingOp(Hk,N,L=L,K=K)
    
    # return
    return OUT
end

function weighting_function(f_cost)

    update_weights = if f_cost == "tukey"

        function tukey_weights(W,r,β)

            # MADNESS (!!!)
            MAD = median( abs.( r .- median( abs.(r) ) ) );

            # Scale
            γ = 1.4825*MAD;

            @inbounds for i in eachindex(r)
                flag = abs(r[i] ./ γ) <= β
                W[i] = flag ? (1 - (abs(r[i])/γ/β)^2 )^2 : zero(r[i]);
            end
        end
        tukey_weights        

    elseif f_cost == "hubber"
        println("No can do, not just yet...")
    elseif f_cost == "google"
        println("No can do, not just yet...")
    elseif f_cost == "lorenz"
        println("No can do, not just yet...")
    end

end

##################################################################
# BFGD-based SSA function for L_2 minimization
function BFGDSSAOp(IN::AbstractArray{T,M},
                   k::Int,
                   maxIter::Int,
                   λ::Real;
                   init::Int = 2) where {T <: Number,M}
# SSA Operator via Bifactored Gradient Descent
# The rank-reduction operation is done via bifactored regression on UV'.
# Syntax:
#       OUT=BFGDSSAOp(IN,PARAM)
# INPUTS:
#       IN - Monochromatic frequency slice
#        k - Desited rank
# OUTPUTS:
#      OUT - Filtered frequency slice
#
# Breno Bahia
##

    # Define sampling
    FwdS,AdjS = Sampling(IN);
    # TODO: Sf,Sa = Sampling(IN)
    # where IN could be Vector, Matrix, Tensor
    # (AbstractArray) but there are two different sampling functions
    # one that removes missing samples and its "put back" adjoint
    # and another that is point-wise multiplication by 1 or 0
    # which should be orthogonal.
    
    # Hankel matrix dimensions
    N = size(IN);
    L = Int.(floor.(N ./ 2) .+ 1);
    K = Int.(N .- L .+ 1);
   
    # Sample the observed data (remove zeroed sample)
    IN0 = FwdS(IN);
    
    # Define forward and adjoint operators
    AOp(x) = FwdS( AveragingOp(x,N,L=L,K=K) ) # Forward operator S(A(X))
    HOp(x) = HankelOp( AdjS(x),N=N,L=L,K=K )  # Adjoint operator A'(S'(X))
    
    # Initialize initial guess
    U0,V0,η = bfgd_init(AOp,HOp,IN0,prod(L),prod(K),k; init=init)

    # Call bfgd for IN0
    Hk,_ = bfgd!(AOp, HOp, IN0, U0, V0, η;
                 λ = λ, maxIter = maxIter,
                 verbose=false, tol=1e-5)
    
    # frequency slice (reconstruct everywhere)
    OUT = AveragingOp(Hk,N,L=L,K=K)
    
    # return
    return OUT
end

##################################################################
# Initialization
function bfgd_init(Fwd, Adj, IN, N, M, k; init::Int = 1)

    # size
    n = (N,M);
    
    # gradient handle
    f_grad(X) = -( Adj( IN .- Fwd( X ) ) );

    # grads at zero and one
    G0 = f_grad( zeros(n) );
    G1 = f_grad(  ones(n) );
    
    # μ estimated (Eq. X in Park et. al. 2018)
    μ = 2*norm(G0-G1,2)/minimum(n);

    # Random initialization
    if init == 0

        # Normalized factors U,V
        Ur = randn(N,k); Ur .= Ur ./ norm(Ur,2);
        Vr = randn(M,k); Vr .= Vr ./ norm(Vr,2);

        # Initial guess and gradient
        Xr = Ur*Vr';
        Gn = f_grad(Xr);        

        # For step-size
        _,σ1,_ = tsvd([Ur; Vr],1)
        _,σ2,_ = tsvd(Gn,1)
        
    elseif init == 1
        
        # Initial guess
        Xr = -(1/μ) * G0;
        
        # Factorize Xr
        Ur,σ1,Vr = tsvd(Xr,k);
        Ur .= Ur * diagm( 0 => sqrt.(σ1) ); Ur .= Ur ./ norm(Ur,2);
        Vr .= Vr * diagm( 0 => sqrt.(σ1) ); Vr .= Vr ./ norm(Vr,2);
        Xr .= Ur * Vr'; 
        Xr .= Adj(Fwd(Xr))
        
                  
        # Step-size
        Gn = f_grad(Xr)
        _,σ2,_ = tsvd(Gn,1)
        
    elseif init == 2 # My initialization
        
        # Hankel matrix from data
        Xr = Adj( IN ) # This is the same as Xr above but not scaled to -1/μ
        
        # tsvd
        Ur,σ1,Vr = tsvd(Xr,k);
        
        # Normalize
        Ur .= Ur * diagm( 0 => sqrt.(σ1) ); Ur .= Ur ./ norm(Ur,2);
        Vr .= Vr * diagm( 0 => sqrt.(σ1) ); Vr .= Vr ./ norm(Vr,2);
        Xr .= Ur * Vr';

        # Make it Hankel
        Xr .= Adj(Fwd(Xr))
        
        # For step-size
        Gn = f_grad(Xr)
        _,σ2,_ = tsvd(Gn,1)
    end
    
    if init == 0
        η = Float32(1 / (0.1 * μ * σ1[1] + σ2[1]));
    else
        η = Float32(1 / (2 * μ * σ1[1] + σ2[1]));
    end

    return Ur,Vr,η
end
##################################################################
