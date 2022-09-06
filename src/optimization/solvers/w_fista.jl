import Base: iterate
import Base.Iterators: take, enumerate

export WFISTAIterable, WFISTAState, wfista_iterable, wfista!

# FISTAIterable
mutable struct WFISTAIterable{Fwd, Adj, Tt, Tw, P, Td, Tx, T <: Number}
    # Operators
    L::Fwd
    Lt::Adj

    # Threshold
    thresh::Tt # This is a vector with a threshold schedule 

    # Update weights
    update_weights::Tw
    args::P

    # Observed data
    d::Td
    
    # Solution
    x0::Tx

    # Things to stop iterating
    α::T
    tol::Number
end

# State
mutable struct WFISTAState{Td, Tx, T <: Number}
    # Arrays
    x::Tx
    xt::Tx
    r::Td
    g::Tx
    W::Td
    
    # Step-size
    t::T
    
    # Iteration related
    tol::Real
    curr_misfit::Real
    prev_misfit::Real
    gprod::Real
    count::Int
end

# Constructor
function wfista_iterable(L, Lt, d, x, thresh, update_weights, args...  ;
                         α = one( real( eltype(d) ) ),
                         tol = sqrt( eps( real( eltype(x) ) ) ),
                         kwargs... )

    return WFISTAIterable(L, Lt, thresh, update_weights, args..., d, x, α, tol)
end

# iterate overload for zero iteration
function iterate(iter::WFISTAIterable{Fwd, Adj, Tt, Tw, P, Td, Tx, T}) where {Fwd, Adj, Tt, Tw, P, Td, Tx, T }

    # count
    count = 1;
    
    # Initial guess
    x = copy(iter.x0)
    xt = copy(x)

    # Residual
    r = iter.d .- iter.L(x)

    # Weights
    W = zero(r);
    iter.update_weights(W,r,iter.args, count)

    # Gradient
    g = iter.Lt(W.*r)

    # Useful quantities
    t = 1.0;
    gprod = real( dot(g,g) )
    prev_misfit = real( dot(iter.d,iter.d) );
    curr_misfit = real( dot(r,r) );

    # define state
    state = WFISTAState{Td, Tx, T}(x, xt, r, g, W, t, iter.tol, gprod, curr_misfit, prev_misfit, count)

    return state,state
end

# subsequent iterations
function iterate(iter::WFISTAIterable{Fwd,Adj,Tt,Tw,P,Td,Tx,T}, state::WFISTAState{Td,Tx,T}) where {Fwd, Adj, Tt, Tw, P, Td, Tx, T}

    # update iteration count
    state.count += 1

    # tmp variables
    tmpx = state.x;
    tmpt = state.t;

    # update on iter.x
    state.x .= state.xt .+ state.g ./ iter.α

    # threshold
    Tα = iter.thresh[state.count]/(2*iter.α)
    soft_threshold!(state.x,Tα)

    # old step-size
    state.t = (1+sqrt(1+4*state.t^2))/2;
    β = (tmpt-1)/state.t;

    # Actual update
    state.xt .= state.x .+ β .* (state.x-tmpx);
    
    # residual
    state.r .= iter.d .- iter.L(state.xt)
    state.prev_misfit = state.curr_misfit
    state.curr_misfit = real( dot(state.r, state.r) )

    # update weights
    iter.update_weights(state.W,state.r,iter.args, state.count)
    
    # gradient
    state.g = iter.Lt(state.W .* state.r);
    state.gprod = dot(state.g, state.g);
   
    # return
    return state,state
end

# These macros will help with a single `iterate` overload
@inline converged(state::WFISTAState) = abs((state.curr_misfit-state.prev_misfit)/state.curr_misfit) <= state.tol

function wfista!(L, Lt, d, x, thresh, update_weights, args... ;
                 α = one( real( eltype(x) ) ),
                 tol = sqrt( eps( real( eltype(x) ) ) ),
                 maxIter::Int = max(1000, length(x) ),
                 verbose::Bool = true,
                 kwargs ... )
    
    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = tol;          
    reserve!(hist, :misfit, maxIter+1)
    
    # define iterable
    iter = wfista_iterable(L,
                           Lt,
                           d,
                           x,
                           thresh,
                           update_weights,
                           args...;
                           α = α,
                           tol = tol,
                           kwargs...)

    iter = halt(iter,converged)
    iter = take(iter,maxIter)

    # This is the loop function
    out = nothing
    for (it, state) in enumerate(iter)
        out = state;
        nextiter!(hist)
        push!(hist, :misfit, out.curr_misfit)
        verbose && println("Iteration $(it) misfit $(out.curr_misfit)")
    end

    shrink!(hist)

    # return
    return out.x, hist
end
