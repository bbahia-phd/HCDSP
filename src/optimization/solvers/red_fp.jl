import Base: iterate
import Base.Iterators: take, enumerate

export REDFPIterable, REDFPState, redfp_iterable, red_fp!

# Iterable
struct REDFPIterable{Op, Rhs, F, Tp, P, Tx, T}
    # Operators
    LOp::Op
    rhs::Rhs
    cost_f::F
    
    # Projection operator
    proj::Tp
    args::P
    
    # Solution & ideal
    x0::Tx
    xi::Tx

    # internal iterations
    int_it::Int

    # tolerance
    ε::T
end

# State
mutable struct REDFPState{Tx, T <: Real}
    # Arrays
    x::Tx     # current solution
    xo::Tx    # aux to store old solution

    # Iteration related
    ε::T     # for δ_new > ε * δ_old
    δ0::T    # zero residual energy
    δ_new::T # curr residual energy
    δ_old::T # prev residual energy
    snr::Any  # for quality
    it::Int  # Iteration counter

end

# Constructor
function redfp_iterable(L, Lt, d, x, μ, proj, args...;
                        ideal = nothing,
                        int_it::Int = max(1000, length(x)),
                        ε = sqrt( eps( real( eltype(x) ) ) ),
                        kwargs...)
    # Define rhs
    aux = Lt(d);
    rhs(δ) = aux .+ μ .* δ;

    # Define Op
    LOp(δ) = Lt(L(δ)) .+ μ .* δ;

    # Define cost function
    cost_f(δ,δ_d) = 0.5 * real(norm(d .- L(δ),2))^2 #+ 0.5 * μ * real(dot(δ,(δ - δ_d)))
    return REDFPIterable(LOp, rhs, cost_f, proj, args..., x, ideal, int_it, ε)
end

# Iterate overload for zero iteration
function iterate(iter::REDFPIterable{Op, Rhs, F, Tp, P, Tx, T}) where {Op,Rhs,F,Tp,P,T,Tx <: AbstractArray}

    # counter
    it = 0;

    # initial solution
    x  = copy(iter.x0)    
    xo = copy(iter.x0)

    δ0    = iter.cost_f(x,x);
    δ_new = copy(δ0)
    δ_old = zero(δ0)
  
    snr = prediction_quality(x,iter.xi);
    
    # define state
    state = REDFPState{Tx,T}(x, xo, iter.ε, δ0, δ_new, δ_old, snr, it)

    return state,state
end

# subsequent iterations
function iterate(iter::REDFPIterable{Op, Rhs, F, Tp ,P, Tx, T}, state::REDFPState{Tx,T}) where {Op,Rhs,F,Tp,P, T, Tx <: AbstractArray}

    # counter
    state.it += 1
    
    # projection
    tmp = iter.proj(state,iter.args)

    # rhs
    tmp .= iter.rhs(tmp);

    # model update
    state.xo .= state.x
    state.x,_ = cg!(iter.LOp,
                    tmp,
                    zero(tmp),
                    max_iter = iter.int_it,
                    verbose=false,
                    tol=T(1e-16))

    # snr (this can be an iterable)
    if iter.xi != nothing;
        state.snr = prediction_quality(state.x,iter.xi);
    end

    # misfit
    state.δ_old = state.δ_new
    state.δ_new = iter.cost_f(state.x,tmp)

    # return 
    return state, state
end

# Convergence
#@inline converged(state::REDFPState) = abs((state.δ_new_misfit-state.δ_old)/state.δ_old) <= state.ε
@inline converged_res(state::REDFPState) = state.δ_new < state.ε^2 * state.δ0
@inline converged_mod(state::REDFPState) = state.it > 0 ? norm(state.xo .- state.x,2)^2/norm(state.xo,2)^2 < state.ε : false;

function red_fp!(L, Lt, d, x, μ, proj, args...;
                 ideal = nothing,
                 ε = sqrt( eps( real( eltype(d) ) ) ),
                 max_iter_o::Int = max(1000, length(x) ),
                 max_iter_i::Int = max(1000, length(x) ),
                 verbose::Bool = true,
                 kwargs ... )

    # Initialize summary
    hist = IterationHistory()
    hist[:tol] = ε
    reserve!(hist, :misfit, max_iter_o+1)
    reserve!(hist, :snr   , max_iter_o+1)
    
    # define iterable
    iter = redfp_iterable(L, Lt, d, x, μ,
                          proj, args...;
                          ideal = ideal,
                          int_it = max_iter_i,
                          ε = ε,
                          kwargs...)
    
#    iter = halt(iter,converged_res)
    iter = halt(iter,converged_mod)    
    iter = take(iter,max_iter_o)

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
