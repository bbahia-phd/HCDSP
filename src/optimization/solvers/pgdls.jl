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
    xo::Tx # previous solution
    x::Tx  # current solution
    z::Tx  # current projected solution
    
    r::Td  # has to match data space
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
    z  = copy(iter.x0)
    x  = copy(iter.x0)
    xo = copy(iter.x0)

    # data and residuals
    r = similar(iter.d);
    
    # model data and residuals
    iter.L(r,x);
    r .= iter.d .- r;
    
    # associated quantities
    δ_new = real(dot(r,r));
    δ_old = zero(δ_new);

    # Gradient
    g = similar(x);
    iter.Lt(r,g);
    γ = real(dot(g,g))

    # Initial Quality
    snr = prediction_quality(x,iter.xi);
    
    # Define state
    state = PGDLSState{Td,Tx,T}(xo, x, z, r, g, iter.ε, snr, γ, δ_new, δ_old, it)

    return state,state
end

# subsequent iterations
function iterate(iter::PGDLSIterable{Fwd,Adj,Tp,P,Td,Tx,T}, state::PGDLSState{Td,Tx,T}) where {Fwd, Adj, Tp,P, T,Td <: AbstractArray, Tx <: AbstractArray}

    # update iteration count
    state.it += 1
    
    # model update
    state.xo .= state.x
    state.z .= state.x .+ iter.α .* state.g

    # projection
    iter.proj(state,iter.args)

    # backtrack
    backtrack(iter,state; proj = true)

    # snr
    if iter.xi !== nothing;
        state.snr = prediction_quality(state.z,iter.xi);
    end

    # residual update
    iter.L(state.r,state.x);
    state.r .= iter.d .- state.r;

    # misfit
    state.δ_old = state.δ_new
    state.δ_new = real(dot(state.r, state.r))

    # gradient
    iter.Lt(state.r,state.g);
    state.γ = real(dot(state.g, state.g))
    
    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::PGDLSState)     = abs((state.δ_new-state.δ_old)/state.δ_old) <= state.ε
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