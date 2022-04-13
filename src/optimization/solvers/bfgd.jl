import Base: iterate
import Base.Iterators: take, enumerate

export BFGDIterable, BFGDState, bfgd_iterable, bfgd!

# BFGDIterable
mutable struct BFGDIterable{Fwd, Adj, Td, Tx, T}
    # Operators
    L::Fwd
    Lt::Adj

    # Observed data
    d::Td

    # Initialization
    U0::Tx
    V0::Tx

    # Step-size
    η::T

    # Regularization
    λ::T
    
    # Things to stop iterating
    tol::Number
end

# State
mutable struct BFGDState{Td, Tx}
    # Residual
    r::Td

    # Initial factors
    Uk::Tx
    Vk::Tx

    # Iterate and copy
    Xk::Tx
    Xp::Tx

    # Gradients
    G::Tx
    Gu::Tx
    Gv::Tx

    # Iteration related
    tol::Number
    curr_misfit::Real
    prev_misfit::Real
    gprod::Real
    count::Int
end

# Constructor
function bfgd_iterable(L, Lt, d, Uk, Vk, η;
                       λ = Float32(1/100),
                       tol = sqrt( eps( real( eltype(x) ) ) ),
                       kwargs... )

    return BFGDIterable(L, Lt, d, Uk, Vk, η, λ, tol)
end

# iterate overload for iterate
function iterate(iter::BFGDIterable{Fwd, Adj, Td, Tx, T}) where {Fwd, Adj, Td, Tx, T}

    # Count
    count = 1;
    
    # Initial guess
    Uk = copy(iter.U0);
    Vk = copy(iter.V0);
    Xk = Uk * Vk';
    Xp = copy(Xk);

    # Residual
    r = copy(iter.d) #iter.d .- iter.L(Xk);

    # Gradients
    G  = copy(Xk); #iter.Lt(r)
    Gu = copy(Uk); #G  * Vk .+ iter.λ .* Uk*(Uk'*Uk - Vk'*Vk)
    Gv = copy(Vk); #G' * Uk .+ iter.λ .* Vk*(Vk'*Vk - Uk'*Uk)
    
    # Useful quantities
    gprod = T( dot(G,G) )
    curr_misfit = T( dot(r,r) );
    prev_misfit = zero(curr_misfit);
    
    # define state
    state = BFGDState{Td,Tx}(r, Uk, Vk, Xk, Xp, G, Gu, Gv, iter.tol, curr_misfit, prev_misfit, gprod, count)
    
    return state,state
end

# subsequent iterations
function iterate(iter::BFGDIterable{Fwd, Adj, Td, Tx, T}, state::BFGDState{Td,Tx}) where {Fwd, Adj, Tw, Td, Tx, T}

    # update iteration count
    state.count += 1

    # residual
    state.r .= iter.d .- iter.L(state.Xp)
    
    # gradient
    state.G  .= iter.Lt(state.r);
    state.Gu .= state.G  * state.Vk .+ iter.λ .* state.Uk*(state.Uk'*state.Uk - state.Vk'*state.Vk)
    state.Gv .= state.G' * state.Uk .+ iter.λ .* state.Vk*(state.Vk'*state.Vk - state.Uk'*state.Uk)
    state.gprod = T( dot(state.G, state.G) )

    # stop cond on grad
    state.prev_misfit = state.curr_misfit
    state.curr_misfit = T( dot(state.r, state.r) )

    # updates
    state.Uk .+= iter.η .* state.Gu;
    state.Vk .+= iter.η .* state.Gv;
    state.Xk  .= state.Uk*state.Vk';
    state.Xp  .= state.Xk;
    
    # return 
    return state, state
end

# These macros will help with a single `iterate` overload
@inline converged(state::BFGDState) = abs((state.curr_misfit-state.prev_misfit)/state.curr_misfit) <= state.tol

function bfgd!(L, Lt, d, Uk, Vk, η;
              λ = Float32(1/4),
              tol = sqrt( eps( real( eltype(d) ) ) ),
              maxIter::Int = max(1000, length(d) ),
              verbose::Bool = true,
              kwargs ... )

    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = tol;          
    reserve!(hist, :misfit, maxIter+1)
    
    # define iterable
    iter = bfgd_iterable(L,
                         Lt,
                         d,
                         Uk,
                         Vk,
                         η;
                         λ = λ,
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
        verbose && println("Iteration $(it) misfit $(out.curr_misfit) gradient $(out.gprod)")
    end

    shrink!(hist)

    # return
    return out.Xk, hist
end
