import Base: iterate
import Base.Iterators: take, enumerate

export REDSDIterable, REDSDState, redsd_iterable, red_sd!

# Iterable
struct REDSDIterable{Fwd, Adj, Tp, P, Td, Tx, T}
    # Operators
    L::Fwd
    Lt::Adj

    # Projection operator
    proj::Tp
    args::P
    
    # Observed data, residual
    d::Td
    
    # Solution & ideal
    x0::Tx
    xi::Tx

    # Step-size
    α::T

    # Trade-off parameter
    μ::T

    # Tolerance
    tol::T
end

# State
mutable struct REDSDState{Td, Tx, T <: Real}
    # Arrays
    x::Tx
    r::Td # has to match data space
    g::Tx # has to match solution space

    # Iteration related
    tol::T
    snr::T
    gprod::T
    curr_misfit::T
    prev_misfit::T
    count::Int
end

# Constructor
function redsd_iterable(L, Lt, d, x, μ, proj, args...;
                        ideal = nothing,
                        α = one( eltype(x) ),
                        tol = sqrt( eps( real( eltype(x) ) ) ),
                        kwargs...)
                        
    return REDSDIterable(L, Lt, proj, args..., d, x, ideal, α, μ, tol)
end


# Iterate overload for zero iteration
function iterate(iter::REDSDIterable{Fwd,Adj,Tp,P,Td,Tx,T}) where {Fwd,Adj,Tp,P,Td,T,Tx <: AbstractArray{T}}

    count = 1;
    x = copy(iter.x0)
    r = iter.d .- iter.L(x)
    g = iter.Lt(r)
    gprod = dot(g,g)
    curr_misfit = dot(r,r);
    prev_misfit = zero(curr_misfit);
    snr = prediction_quality(x,iter.xi);
    
    # define state
    state = REDSDState{Td,Tx,T}(x, r, g, iter.tol, snr, gprod, curr_misfit, prev_misfit, count)

    return state,state
end

# subsequent iterations
function iterate(iter::REDSDIterable{Fwd,Adj,Tp,P,Td,Tx,T}, state::REDSDState{Td,Tx,T}) where {Fwd, Adj, Tp,P, T,Td <: AbstractArray{T}, Tx <: AbstractArray{T}}

    # counter
    state.count += 1
    
    # projection
    tmp = iter.proj(state,iter.args)

    # gradient
    state.g .= iter.Lt(state.r) .- iter.μ .* (state.x .- tmp);
    state.gprod = dot(state.g, state.g)

    # model update
    #state.x .+= iter.α .* state.g

    # backtrack
    backtrack(iter,state)

    # snr (this can be an iterable)
    if iter.xi != nothing;
        state.snr = prediction_quality(state.x,iter.xi);
    end

    # residual update
    state.r .= iter.d .- iter.L(state.x)
    
    # misfit
    state.prev_misfit = state.curr_misfit
    state.curr_misfit = dot(state.r, state.r)

    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::REDSDState) = abs((state.curr_misfit-state.prev_misfit)/state.curr_misfit) <= state.tol

function red_sd!(L, Lt, d, x, μ, proj, args...;
                ideal = nothing,
                α = one( eltype(x) ),
                tol = sqrt( eps( real( eltype(d) ) ) ),
                max_iter::Int = max(1000, length(x) ),
                verbose::Bool = true,
                kwargs ... )

    # Initialize summary
    hist = IterationHistory()
    hist[:tol] = tol
    reserve!(hist, :misfit, max_iter+1)
    reserve!(hist, :snr, max_iter+1)
    
    # define iterable
    iter = redsd_iterable(L,
                          Lt,
                          d,
                          x,
                          μ,
                          proj,
                          args...;
                          ideal = ideal,
                          α = α,
                          tol = tol,
                          kwargs...)
    
    iter = halt(iter,converged)
    iter = take(iter,max_iter)

    # This is the loop function
    out = nothing
    for (it, state) in enumerate(iter)
        out = state;
        nextiter!(hist)
        push!(hist, :misfit, out.curr_misfit)
        push!(hist, :snr, out.snr)
        verbose && println("Iteration $(it) misfit $(out.curr_misfit) snr $(out.snr)")
    end

    shrink!(hist)

    # return
    return out.x, hist
end
