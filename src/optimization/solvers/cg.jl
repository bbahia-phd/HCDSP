import Base: iterate
import Base.Iterators: enumerate

export CGIterable, CGState, cg_iterable, cg!

# CGIterable
struct CGIterable{TOp, Td, Tx, T <: Real}
    # Operators
    L::TOp
    # Observed data
    d::Td
    # Solution
    x0::Tx
    # Tolerance
    ε::T
end

# State
mutable struct CGState{Ta,Tb}
    # Arrays
    x::Ta     # current solution
    xo::Ta    # aux to store old solution
    r::Ta     # residual
    p::Ta     # conjugate vector
    Ap::Ta    # aux to conjugate vector

    δ0::Tb    # zero residual energy
    δ_new::Tb # curr residual energy
    δ_old::Tb # prev residual energy
    ε::Tb     # for δ_new > ε * δ_old

    it::Int   # Iteration counter
end

# Constructor
function cg_iterable(L, d, x ;
                     ε = sqrt( eps( real( eltype(x) ) ) ))
    return CGIterable(L, d, x, ε)
end


# iterate overload
function iterate(iter::CGIterable{TOp, Td, Tx, T}) where {TOp, Td, T, Tx}
    
    # initial solution
    x  = copy(iter.x0)    
    xo = copy(iter.x0)
    
    # residual
    r = iter.d - iter.L(xo)

    # Conjugate vecors
    p = copy(r);
    Ap = similar(r);

    # To compute step-sizes
    δ0 = real(dot(r,r));
    δ_new = copy(δ0);
    δ_old = zero(δ0);
    
    # Iteration counter
    it = 0;

    # Call state
    state = CGState{Tx,T}(x, xo, r, p, Ap, δ0, δ_new, δ_old, iter.ε, it)

    # return 
    return state, state
end

function iterate(iter::CGIterable{TOp, Td, Tx, T}, state::CGState{Ta, Tb}) where {TOp, Td, T, Tx <: AbstractArray, Ta, Tb}

    # Update
    state.it += 1;
    
    # Single activation
    state.Ap .= iter.L(state.p);

    # Step-size
    α = state.δ_new / real(dot(state.p, state.Ap))

    # Updates
    state.xo .= state.x
    state.x .+= α .* state.p
    state.r .-= α .* state.Ap

    # Step-size
    state.δ_old = state.δ_new
    state.δ_new = real(dot(state.r, state.r))
    β = state.δ_new / state.δ_old

    # Update
    state.p .= state.r .+ β .* state.p

    # Return
    return state, state
end

# Convergence
@inline converged_res(state::CGState) = state.δ_new < state.ε^2 * state.δ0
@inline converged_mod(state::CGState) = state.it > 0 ? norm(state.xo .- state.x,2)^2/norm(state.xo,2)^2 < state.ε : false;

function cg!(L, d, x ;
             tol = sqrt( eps( real( eltype(d) ) ) ),
             max_iter::Int = max(1000, length(x) ),
             verbose::Bool = true,
             kwargs ... )

    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = tol;
    reserve!(hist, :misfit, max_iter+1)
    
    # define iterable
    iter = cg_iterable(L, d, x; ε = tol, kwargs...)
    iter = halt(iter, converged_res)
    iter = halt(iter, converged_mod)    
    iter = take(iter, max_iter)

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
