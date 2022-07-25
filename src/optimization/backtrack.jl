import Base: iterate
import Base.Iterators: take, enumerate
import NaNMath: min,max
    
export backtrack, backtrack_iterable, BTIterable, BTState, get_ϕ, PBTIteable, PBTState


#####
##### Backtracking iterable
#####
"""
The `backtrack` function receives an iterable and is associated state.
It is assumed that both follow the solvers coded in this package to
solve least-squares problems.

Therefore, it looks for specific fields such as the forward operator,
the gradient energy and the cost function at a given iteration.

In addition, it defaults some constants used in inexact line searches,
here in specific the ones used in Armijo's rule which are usually denoted by α_0, σ_1, γ_1, γ_2.


May 26, 22: I added real(dot(x,y)) in some cases to help backtrack in the case of quaternion variables.
Despite the quaternion numbers, the step-sizes are obtained as real numbers (similarly to the complex case).

"""
# calling function
function backtrack(iter, state; 
                   α_0  = one( eltype(state.γ) ),
                   σ_1  = convert.( typeof(α_0), 1e-4 ),
                   proj = nothing,
                     W  = nothing,
                   max_iter = 1000,
                   verbose::Bool=false)

    # Allocate output
    out = zero(state.x)

    # Initialize backtrack
    bt_iter = backtrack_iterable(iter,
                                 state,
                                 out;
                                 W = W,
                                 proj = proj,
                                 α_0 = α_0,
                                 σ_1 = σ_1,
                                 max_iter = max_iter);
    bt_iter = halt(bt_iter, sufficient_decrease);
    bt_iter = take(bt_iter, max_iter)

    # Backtrack iterations
    aux = nothing
    for (it, bt_state) in enumerate(bt_iter)
        aux = bt_state;
        verbose && println("Backtrack step-size $(aux.α_c) at iteration $(it).")
    end

    # update solution
    if proj == nothing
        state.x .+= aux.α_c .* state.g;
    else
        state.x .+= aux.α_c .* (proj-state.x);
    end
end

# constructor
function backtrack_iterable(iter, state, out;
                            W = nothing,
                            proj = nothing,
                            α_0 = one( real( eltype( state.x  ) ) ),
                            σ_1 = convert.( typeof(α_0), 1e-4 ),
                            γ_1 = convert.( typeof(α_0), 0.1 ),
                            γ_2 = convert.( typeof(α_0), 0.5 ),
                            max_iter = 1000,
                            order=3)

    if proj == nothing
        # Function handle for LS misfit as function of step-size α
        ϕ = get_ϕ(iter.L, iter.d, state.x, state.g, out; W = W)
        return BTIterable(ϕ,
                          state.γ,
                          state.δ_new,
                          α_0,
                          σ_1,
                          γ_1,
                          γ_2,
                          max_iter,
                          order)
    else
        # Function handle for LS misfit as function of step-size α
        ϕ = get_ϕ(iter.L, iter.d, state.x, (proj-state.x), out; W = W)
        return BTIterable(ϕ,
                           real(dot(state.g,(proj-state.x))),
                           state.δ_new,
                           α_0,
                           σ_1,
                           γ_1,
                           γ_2,
                           max_iter,
                           order)
    end
end

"""
A function handle for a quadratic cost function
"""
function get_ϕ( L, d, x, g, xnew ;
                r = zero(d),
                W = nothing )
    ϕ = if W == nothing
        function _ϕ(α)
            # Step towards g
            xnew .= x .+ α .* g;

            # Residual
            r .= d .- L(xnew);

            # Sum of squares
            return real(dot(r,r))
        end
        _ϕ
    else
        function _ϕw(α)
            # Step towards g
            xnew .= x .+ α .* g;

            # Residual
            r .= W .* (d .- L(xnew));

            # Sum of squares
            return real(dot(r,r))
        end
        _ϕw
    end
    return ϕ
end

struct BTIterable{F,T}
    # Function handle
    ϕ::F
    
    # Constants
    dϕ_0::T # gprod
    ϕ_0::T  # current misfit
    α_0::T  # initial step-size
    σ_1::T  # reduction constant
    γ_1::T  # lower safeguard
    γ_2::T  # upper safeguard
    max_iter::Int # number of iterations
    order::Int    # order of interpolants
end

mutable struct BTState{T}
    # candidate step-size
    α_c::T

    # previous step-size
    α_p::T

    # misfit and slope
    ϕx0::T # misfit for α_p (used for cubic interpolation)
    ϕx::T  # misfit for α_c
    dϕ::T  # sufficient decrease

    # counter
    count::Int
end

# what is in bt_state
function iterate(iter::BTIterable{F,T}) where {F,T}

    # This variable holds the number of max iterations
    # the code should perform to find a finite value
    # for the function being minimized.
    finite_max_iter = -log2( eps( real( T ) ) )
    
    # initialize
    ϕx_1 = iter.ϕ_0
    α_p, α_c  = iter.α_0, iter.α_0

    # Evaluate cost function using initial step-size
    ϕx_1 = iter.ϕ(α_p)

    # Trivial backtrack to find finite function value
    count_finite = 0;
    while !isfinite(ϕx_1) && count_finite < finite_max_iter
        # count
        count_finite += 1

        # store and halve
        α_p = α_c
        α_c = α_p .* 0.5

        # eval cost func with halved step-size
        ϕx_1 = iter.ϕ(α_c)
    end

    # slope for sufficient decrease condition
    dϕ = iter.ϕ_0 + α_c .* iter.σ_1 * iter.dϕ_0

    # define state
    state = BTState{T}(α_c, α_p, ϕx_1, ϕx_1, dϕ, 0)

    # return
    return state,state      
end

# what is in bt_state
function iterate(iter::BTIterable{F,T}, state::BTState{T}) where {F,T}

    # count 
    state.count += 1
    
    # backtrack step-size
    if iter.order == 2 || state.count == 1
        # quadratic interpolant
        α_tmp = -(iter.dϕ_0 * state.α_c^2)/(2*(state.ϕx - iter.ϕ_0 - iter.dϕ_0 * state.α_c))
    else
        # cubic interpolant
        div = one(T) / (state.α_p^2 * state.α_c^2 * (state.α_c - state.α_p))
        a = ( state.α_p^2*(state.ϕx - iter.ϕ_0 - iter.dϕ_0*state.α_c ) - state.α_c^2*( state.ϕx0 - iter.ϕ_0 - iter.dϕ_0*state.α_p  )) * div
        b = (-state.α_p^3*(state.ϕx - iter.ϕ_0 - iter.dϕ_0*state.α_c ) + state.α_c^3*( state.ϕx0 - iter.ϕ_0 - iter.dϕ_0*state.α_p  )) * div

        if isapprox(a, zero(a), atol=eps(real(T)))
            α_tmp = iter.dϕ_0 / (2*b);
        else
            # delta
            d = max(b^2 - 3*a*iter.dϕ_0, T(0))
            # zeros
            α_tmp = (-b + sqrt(d)) / (3*a)
        end
    end

    # update previous step-size
    state.α_p = state.α_c

    # Avoid too small/big corrections (see package NaNMath.jl)
        α_tmp = NaNMath.min(α_tmp, state.α_c * iter.γ_2)
    state.α_c = NaNMath.max(α_tmp, state.α_c * iter.γ_1)

    # candidate misfit
    state.ϕx0 = state.ϕx
    state.ϕx  = iter.ϕ(state.α_c)

    # slope for sufficient decrease condition
    state.dϕ  = iter.ϕ_0 + state.α_c .* iter.σ_1 * iter.dϕ_0
   
    return state,state      
end

# Armijo condition
@inline sufficient_decrease(state::BTState) = real(state.ϕx) <= real(state.dϕ)
