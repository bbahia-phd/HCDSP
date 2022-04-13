import Base: iterate
import Base.Iterators: take, enumerate

export IHTIterable, IHTState, iht_iterable, iht!, threshold!

# IHTIterable
mutable struct IHTIterable{Fwd, Adj, Tt, Td, Tx, T <: Number}
    # Operators
    L::Fwd
    Lt::Adj

    # Threshold
    thresh::Tt # This is a vector with a threshold schedule 

    # Observed data
    d::Td
    
    # Solution
    x0::Tx

    # Things to stop iterating
    alpha::T
    tol::Number
end

# State
mutable struct IHTState{Td, Tx, T <: Number}
    # Arrays
    x::Tx
    r::Td
    g::Tx

    # Iteration related
    tol::Real
    curr_misfit::Real
    prev_misfit::Real
    gprod::Real
    count::Int
end

# Constructor
function iht_iterable(L, Lt, d, x, thresh ;
                      alpha = one( real( eltype(d) ) ),
                      tol = sqrt( eps( real( eltype(x) ) ) ),
                      kwargs... )

    return IHTIterable(L, Lt, thresh, d, x, alpha, tol)
end

# iterate overload for zero iteration
function iterate(iter::IHTIterable{Fwd, Adj, Tt, Td, Tx, T}) where {Fwd, Adj, Tt, T, Td, Tx }

    # Initial guess
    x = copy(iter.x0)
    threshold!(x,iter.thresh[1])

    # Residual and gradient
    r = iter.d .- iter.L(x)
    g = iter.Lt(r)

    # Useful quantities
    count = 0;
    gprod = real( dot(g,g) )
    prev_misfit = real( dot(iter.d,iter.d) );
    curr_misfit = real( dot(r,r) );

    # define state
    state = IHTState{Td, Tx, T}(x, r, g, iter.tol, gprod, curr_misfit, prev_misfit, count)

    return state,state
end

# subsequent iterations
function iterate(iter::IHTIterable{Fwd,Adj,Tt,Td,Tx,T}, state::IHTState{Td,Tx,T}) where {Fwd, Adj, Tt,T,Td, Tx}

    # update iteration count
    state.count += 1

    # auxiliary variable
    q = iter.L(state.g);

    # step-size
    α = T( state.gprod / (dot(q,q) + 1e-13) )    

    # backtrack
    backtrack(iter, state; α_0 = α)

    # threshold
    threshold!(state.x, iter.thresh[state.count])

    # residual
    state.r .= iter.d .- iter.L(state.x)
    state.prev_misfit = state.curr_misfit
    state.curr_misfit = real( dot(state.r, state.r) )

    # gradient
    state.g = iter.Lt(state.r);
    state.gprod = dot(state.g, state.g)
   
    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::IHTState) = abs((state.curr_misfit-state.prev_misfit)/state.curr_misfit) <= state.tol

function iht!(L, Lt, d, x, thresh ;
              alpha = one( real( eltype(x) ) ),
              tol = sqrt( eps( real( eltype(x) ) ) ),
              maxIter::Int = max(1000, length(x) ),
              verbose::Bool = true,
              kwargs ... )

    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = tol;          
    reserve!(hist, :misfit, maxIter+1)
    
    # define iterable
    iter = iht_iterable(L, Lt, d, x, thresh; alpha = alpha, tol = tol, kwargs...)
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

function threshold!(x::AbstractArray{<:Number}, t)
# hard-threshold
    @assert t >= 0
    @inbounds for i in eachindex(x)
        if abs(x[i]) <= t
            x[i] = zero(x[i])
         end
     end
end
