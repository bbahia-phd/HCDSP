import Base: iterate
import Base.Iterators: enumerate

export CGLSIterable, CGLSState, cgls_iterable, cgls!

# CGLSIterable
mutable struct CGLSIterable{TOp, Td, Tx, T <: Real}
    # Operators
    L::TOp
    Lt::TOp

    # Observed data, residual
    d::Td
    r::Td
    
    # Solution, gradient, conj. grad
    x::Tx
    g::Tx
    p::Tx

    # Things to stop iterating
    tol::T
    gprod::T
    misfit::T
    maxIter::Int
end

# State
mutable struct CGLSState{Ta}
    # Arrays
    r::Ta
    g::Ta
    p::Ta
end

# constructor
function cgls_iterable(L, Lt, d, x ;
                       tol = sqrt( eps( real( eltype(x) ) ) ),
                       maxIter::Int = max(1000, length(x) ),
                       state::CGLSState = CGLSState( similar(d), zero(x), zero(x) ))

    # Extract state variables
    r = state.r;
    g = state.g;
    p = state.p;
    
    # Initialize arrays (model, residual, gradient, conj gradient)
    r = d .- L*x;
    g = Lt*r;
    p = copy(g);

    # Scalars (step-sizes, tolerances, and norms)
    gprod = dot(g,g);
    misfit = dot(r,r) * tol ;

    return CGLSIterable(L, Lt, d, r, x, g, p, tol, gprod, misfit, maxIter)
end

# These macros will help with a single `iterate` overload
@inline converged(iter_obj::CGLSIterable) = iter_obj.misfit <= iter_obj.tol
@inline first(iter_obj::CGLSIterable) = 0;
@inline stop(iter_obj::CGLSIterable, iter::Int) = iter >= iter_obj.maxIter || converged(iter_obj);

# iterate overload
function iterate(iter::CGLSIterable{TOp,Td,Tx,T}, it::Int=first(iter)) where {TOp, Td, T, Tx <: AbstractArray{T}}
    
    if stop(iter, it)
        return nothing
    end

    # aux variables
    q = similar(iter.d)
    mul!(q, iter.L, iter.p)

    # alpha
    alpha = iter.gprod / (dot(q,q) + 1e-10)

    # update
    iter.x .+= alpha .* iter.p
    iter.r .-= alpha .* q
    iter.misfit = dot(iter.r,iter.r)

    # gradient
    mul!(iter.g, iter.Lt, iter.r)

    # beta
    gprod_new = dot(iter.g, iter.g)
    beta = gprod_new / (iter.gprod + 1e-10)
    iter.gprod = gprod_new
    
    # conjugate gradient
    iter.p .= iter.g .+ beta .* iter.p

    # return 
    return iter, it + 1
end

################# cgls function ####################
# TODO: add a method to start zero-valued models
# cgls(L,Lt,d; kwargs...) =
# cgls(L,Lt,d,zeros(???); initial_zero = true,  kwargs ...)
# zeros(Lt * d) is a solution: it  should have the size of Lt*d
# and its corresponding type. For now, we set it up from the
# calling environment.
####################################################

function cgls!(L, Lt, d, x ;
               tol = sqrt( eps( real( eltype(d) ) ) ),
               maxIter::Int = max(1000, length(x) ),
               state::CGLSState = CGLSState( similar(d), zero(x), zero(x) ),
               verbose::Bool = true,
               kwargs ... )

    # Initialize summary
    hist = IterationHistory(); # Initialize history
    hist[:tol] = tol;          # Start filling in its fields (notice
    # that this goes within hist.data[:tol], but you can define it
    # like that because we have added methods to do that (see common_iterables/iterables.jl).
    reserve!(hist, :misfit, maxIter+1)
    
    # define iterable
    iter = cgls_iterable(L, Lt, d, x; tol = tol, maxIter = maxIter, state = state, kwargs...)

    for (it, out) in enumerate(iter)
        nextiter!(hist)
        push!(hist, :misfit, iter.misfit)
        verbose && println("Iteration $(it) misfit $(iter.misfit)")
    end

    shrink!(hist)

    # return
    return iter.x, hist

    # outputs should be something like
    # state.x and history
    # history is a tuple with everything you
    # should need from the inversion process
    # check IterativeSolvers.jl
    # to have a better idea
end
