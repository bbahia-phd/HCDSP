import Base: iterate
import Base.Iterators: take, enumerate

export REDADMMIterable, REDADMMState, red_admm_iterable, red_admm!

# Iterable
struct REDADMMIterable{Cf, Gf, Op, Tp, P, Tx, T}
    # Function handles
    cost_f::Cf
    get_grad::Gf
    Op::Op
    
    # Denoiser
    proj::Tp
    args::P

    # Solution & ideal
    x0::Tx
    xi::Tx

    # internal iterations
    it1::Int
    it2::Int

    # Hyperparameters
    μ::T
    β::T
    
    # Tolerance
    tol::T
end

# State
mutable struct REDADMMState{Tx, T <: Real}
    # Arrays
    x::Tx
    v::Tx
    u::Tx
    z::Tx
    zs::Tx

    # Iteration related
    tol::T
    snr::Any
    curr_misfit::T
    prev_misfit::T
    it::Int
end

# Constructor
function red_admm_iterable(L, Lt, d, x, μ, β, proj, args...;
                           ideal = nothing,
                           it1::Int = max(1000, length(x)),
                           it2::Int = max(1000, length(x)),                           
                           tol = sqrt( eps( real( eltype(x) ) ) ),
                           kwargs...)

    # Define cost function
    cost_f(δ,δ_d) = 0.5 * real(norm(d .- L(δ),2))^2 + 0.5 * μ * real( dot(δ,(δ - δ_d)) )

    # Define get_grad
    get_grad(δ,δ_d) = Lt(L(δ)-d) .+ β .* (δ .- δ_d)
    
    # Define Op
    Op(δ) = Lt(L(δ)) .+ β .* δ;

    return REDADMMIterable(cost_f, get_grad, Op, proj, args..., x, ideal, it1, it2, μ, β, tol)
end

# Iterate overload for zero iteration
function iterate(iter::REDADMMIterable{Cf, Gf, Op, Tp, P, Tx, T}) where {Cf, Gf, Op, Tp, P, T, Tx <: AbstractArray}

    it = 1;
    x  = copy(iter.x0)
    v  = copy(iter.x0)
    u  = zero(iter.x0)
    z  = zero(iter.x0)
    zs = zero(iter.x0)

    curr_misfit = iter.cost_f(x,x);
    prev_misfit = zero(curr_misfit);
    snr = prediction_quality(x,iter.xi);
    
    # define state
    state = REDADMMState{Tx,T}(x, v, u, z, zs, iter.tol, snr, curr_misfit, prev_misfit, it)

    return state,state
end

# subsequent iterations
function iterate(iter::REDADMMIterable{Cf, Gf, Op, Tp, P, Tx, T}, state::REDADMMState{Tx, T}) where {Cf, Gf, Op, Tp, P, T, Tx <: AbstractArray}

    # counter
    state.it += 1

    # allocate temps
    e = zero(state.x);
    g = zero(state.x);

    state.z  .= state.x;
    state.zs .= state.v - state.u;

    # repelace by CG call 
    for i in 1:iter.it1
        g .= iter.get_grad(state.z,state.zs);
        e .= iter.Op(g);
        α = real( dot(g,g) / dot(g,e) )
        state.z .-= α .* g;
    end
    
    state.x  .= state.z;
    state.z  .= state.v; 
    state.zs .= state.x .+ state.u;

    for i in 1:iter.it2
        e .= iter.proj(state,iter.args)
        state.z .= (iter.μ .* e .+ iter.β .* state.zs) ./ (iter.β + iter.μ);
    end

    state.v .= state.z;
    state.u .+= (state.x .- state.v);
    
    # snr (this can be an iterable)
    if iter.xi != nothing;
        state.snr = prediction_quality(state.x,iter.xi);
    end

    # misfit
    state.prev_misfit = state.curr_misfit
    state.curr_misfit = iter.cost_f(state.x,state.z)

    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::REDADMMState) = abs((state.curr_misfit-state.prev_misfit)/state.curr_misfit) <= state.tol

function red_admm!(L, Lt, d, x, μ, β, proj, args...;
                   ideal = nothing,
                   tol = sqrt( eps( real( eltype(d) ) ) ),
                   max_iter_o::Int = max(1000, length(x) ),
                   max_iter_i1::Int = max(1000, length(x) ),
                   max_iter_i2::Int = max(1000, length(x) ),
                   verbose::Bool = true,
                   kwargs ... )

    # Initialize summary
    hist = IterationHistory()
    hist[:tol] = tol
    reserve!(hist, :misfit, max_iter_o+1)
    reserve!(hist, :snr   , max_iter_o+1)
    
    # define iterable
    iter = red_admm_iterable(L, Lt, d, x, μ, β,
                             proj, args...;
                             ideal = ideal,
                             it1 = max_iter_i1,
                             it2 = max_iter_i2,
                             tol = tol,
                             kwargs...)
    
    iter = halt(iter,converged)
    iter = take(iter,max_iter_o)

    # This is the loop function
    out = nothing
    for (it, state) in enumerate(iter)
        out = state;
        nextiter!(hist)
        push!(hist, :misfit, out.curr_misfit)
        push!(hist, :snr   , out.snr)
        verbose && println("Iteration $(it) misfit $(out.curr_misfit) snr $(out.snr)")
    end

    shrink!(hist)

    # return
    return out.x, hist
end
