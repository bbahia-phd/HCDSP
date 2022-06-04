import Base: iterate
import Base.Iterators: take, enumerate

export PGDLSIterable, PGDLSState, pgdls_iterable, pgdls!

# CGLSIterable
struct PGDLSIterable{Fwd, Adj, Tp, P, Td, Tx, T}
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

    # step size
    α::T
    
    # tolerance
    ε::T
end

# State
mutable struct PGDLSState{Td, Tx, T <: Number}
    # Arrays
    x::Tx  # current solution
    xo::Tx # old solution
    
    r::Td # has to match data space
    g::Tx # has to match solution space

    # Iteration related
    ε::T   # tol::Number
    snr::Any
    γ::T    #gprod::T
    δ_new::T    #curr_misfit::T
    δ_old::T    #prev_misfit::T
    it::Int    #count::Int
end

# constructor
function pgdls_iterable(L, Lt, d, x, proj, args...;
                        ideal = nothing,
                        α = one( eltype(x) ),
                        ε = sqrt( eps( real( eltype(x) ) ) ),
                        kwargs...)
                        
    return PGDLSIterable(L, Lt, proj, args..., d, x, ideal, α, ε)
end

# iterate overload for zero iteration
function iterate(iter::PGDLSIterable{Fwd,Adj,Tp,P,Td,Tx,T}) where {Fwd,Adj,Tp,P,Td,T,Tx <: AbstractArray}

    # counter
    it = 0;

    # initial solution
    x  = copy(iter.x0)
    xo = copy(iter.x0)

    # residuals
    r = iter.d .- iter.L(x)
    δ_new = dot(r,r);
    δ_old = zero(δ_new);

    # Gradient
    g = iter.Lt(r)
    γ = dot(g,g)

    # Initial Quality
    snr = prediction_quality(x,iter.xi);
    
    # Define state
    state = PGDLSState{Td,Tx,T}(x, xo, r, g, iter.ε, snr, γ, δ_new, δ_old, it)

    return state,state
end

# subsequent iterations
function iterate(iter::PGDLSIterable{Fwd,Adj,Tp,P,Td,Tx,T}, state::PGDLSState{Td,Tx,T}) where {Fwd, Adj, Tp,P, T,Td <: AbstractArray, Tx <: AbstractArray}

    # update iteration count
    state.it += 1
    
    # model update
    state.xo .= state.x
    state.x .+= iter.α .* state.g

    # projection
    tmp = iter.proj(state,iter.args)

    # backtrack
    backtrack(iter,state; proj = tmp)

    # snr
    if iter.xi != nothing;
        state.snr = prediction_quality(state.x,iter.xi);
    end

    # residual update
    state.r .= iter.d .- iter.L(state.x)

    # misfit
    state.δ_old = state.δ_new
    state.δ_new = real(dot(state.r, state.r))

    # gradient
    state.g = iter.Lt(state.r);
    state.γ = real(dot(state.g, state.g))
    
    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::PGDLSState)     = abs((state.δ_new-state.δ_old)/state.δ_old) <=    state.ε
@inline converged_res(state::PGDLSState) = state.δ_new < state.ε^2 * state.δ0
@inline converged_mod(state::PGDLSState) = state.it > 0 ? norm(state.xo .- state.x,2)^2/norm(state.xo,2)^2 < state.ε : false;

function pgdls!(L, Lt, d, x, proj, args...;
                ideal = nothing,
                α = one( eltype(x) ),
                ε = sqrt( eps( real( eltype(d) ) ) ),
                maxIter::Int = max(1000, length(x) ),
                verbose::Bool = true,
                kwargs ... )

    # Initialize summary
    hist = IterationHistory()
    hist[:tol] = ε
    reserve!(hist, :misfit, maxIter+1)
    reserve!(hist, :snr   , maxIter+1)
    
    # define iterable
    iter = pgdls_iterable(L, Lt, d, x,
                          proj, args...;
                          ideal = ideal,
                          α = α,
                          ε = ε,
                          kwargs...)
    iter = halt(iter,converged)
    #    iter = halt(iter,converged_res)
    iter = halt(iter,converged_mod)
    iter = take(iter,maxIter)

    # This is the loop function
    out = nothing
    for (it, state) in enumerate(iter)
        out = state;
        nextiter!(hist)
        push!(hist, :misfit, out.δ_new)
        push!(hist, :snr   , out.snr)
        verbose && println("Iteration $(it) misfit $(out.δ_new) snr $(out.snr)")
    end

    shrink!(hist)

    # return
    return out.x, hist
end