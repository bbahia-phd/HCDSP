import Base: iterate
import Base.Iterators: take, enumerate

export RIHTIterable, RIHTState, riht_iterable, riht!

# IHTIterable
mutable struct RIHTIterable{Fwd, Adj, Tt, Tw, P, Td, Tx, T <: Real}
    # Operators
    L::Fwd
    Lt::Adj

    # Threshold
    thresh::Tt 

    # Update weights
    update_weights::Tw
    args::P

    # Observed data
    d::Td
    
    # Solution & ideal
    x0::Tx

    # Things to stop iterating
    α::T
    ε::Number
end

# State
mutable struct RIHTState{Td, Tx}
    # Arrays
    x::Tx
    r::Td
    g::Tx
    W::Td

    
    # Iteration related
    ε::Number   #tol::Number
    δ_new::Real #curr_misfit::Real
    δ_old::Real #prev_misfit::Real
    γ::Real     #gprod::Real
    it::Int     #count::Int
end

# Constructor
function riht_iterable(L, Lt, d, x, thresh, update_weights, args... ;
                      α = one( real( eltype(d) ) ),
                      ε = sqrt( eps( real( eltype(x) ) ) ),
                      kwargs... )

    return RIHTIterable(L, Lt, thresh, update_weights, args..., d, x, α, ε)
end

# iterate overload for iterate
function iterate(iter::RIHTIterable{Fwd, Adj, Tt, Tw, P, Td, Tx, T}) where {Fwd, Adj, Tt, Tw, P, Td, Tx, T}

    # Count
    it = 1;

    # Initial guess
    x = copy(iter.x0)
    threshold!(x,iter.thresh[it])

    # Residual and gradient
    r = iter.d .- iter.L(x);

    # Weights
    W = zero(r);
    iter.update_weights(W,r,iter.args, it)

    # Gradient
    g = iter.Lt(W .* r)
    γ = T( dot(g,g) )
    
    # Useful quantities
    δ_old = T( dot(r,r) );
    δ_new = copy(δ_old);
    
    # define state
    state = RIHTState{Td,Tx}(x, r, g, W, iter.ε, γ, δ_new, δ_old, it)

    return state,state
end

# subsequent iterations
function iterate(iter::RIHTIterable{Fwd, Adj, Tt, Tw, P, Td, Tx, T}, state::RIHTState{Td,Tx}) where {Fwd, Adj, Tt, Tw, P, Td, Tx, T}

    # update iteration count
    state.it += 1

    # auxiliary variable
    q = sqrt.(state.W) .* iter.L(state.g);

    # step-size
    α = T( state.γ / (dot(q,q) + 1e-13) )    

    # backtrack
    backtrack(iter, state; α_0 = α, W = state.W)

    # threshold
    threshold!(state.x, iter.thresh[state.it])
    
    # residual
    state.r .= iter.d .- iter.L(state.x)
    state.δ_old = T( state.δ_new )
    state.δ_new = T( dot(state.r, state.r) )

    # update weights
    iter.update_weights(state.W,state.r,iter.args,state.it)
    
    # gradient
    state.g .= iter.Lt(state.W .* state.r);
    state.γ = T( dot(state.g, state.g) )
   
    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::RIHTState) = abs((state.δ_new-state.δ_old)/state.δ_new) <= state.ε

function riht!(L, Lt, d, x, thresh, update_weights, args... ;
              α = one( real( eltype(x) ) ),
              ε = sqrt( eps( real( eltype(x) ) ) ),
              maxIter::Int = max(1000, length(x) ),
              verbose::Bool = true,
              kwargs ... )

    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = ε;          
    reserve!(hist, :misfit, maxIter+1)
    
    # define iterable
    iter = riht_iterable(L,
                         Lt,
                         d,
                         x,
                         thresh,
                         update_weights,
                         args...;
                         α = α,
                         ε = ε,
                         kwargs...)
    
    iter = halt(iter,converged)
    iter = take(iter,maxIter)

    # This is the loop function
    out = nothing
    for (it, state) in enumerate(iter)
        out = state;
        nextiter!(hist)
        push!(hist, :misfit, out.δ_new)
        verbose && println("Iteration $(it) misfit $(out.δ_new)")
    end

    shrink!(hist)

    # return
    return out.x, hist
end
