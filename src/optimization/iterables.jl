import Base: iterate, getindex, setindex!, keys
export IterationHistory
    
#=
#= Iterable functionalities to be used by different solvers =#

#####
##### Loop 
#####

# Loop within the iterative method itself

#function loop(iter)
#
#    x = nothing
#    for y in iter
#        x = y;
#    end
#
#    return x
#end
=#

#####
##### Iterable for stop conditions
#####

mutable struct HaltIterable{I,F}
    iter::I # any iterable
    stop_cond::F # a stop condition as a function
end

function iterate(iter::HaltIterable)
    next = iterate(iter.iter)     # The input has a iter field
    return dispatch(iter,next)    # Testing the stop condition
end

function iterate(iter::HaltIterable, (flag, state))
    # check condition
    if flag == :stop
        return nothing
    end

    # Run the iteration using iterable within iterable
    next = iterate(iter.iter,state)
    return dispatch(iter,next)
end

function dispatch(iter::HaltIterable, next)
    # check if iteration is over
    if next === nothing
        return nothing
    end

    # Apply stop_cond to next[1] (actual iterate)
    # and check stop condition
    iter.stop_cond(next[1])
    return next[1], (iter.stop_cond(next[1]) ? :stop : :run, next[2])
end

# constructor
halt(iter::I, stop_cond::F) where {I,F} = HaltIterable{I,F}(iter,stop_cond)

#####
##### Show Iterable
#####

"""
Display summary of methods state.
"""

mutable struct SummaryIterable{I,F}
    iter::I
    summ::F
end

function iterate(iter::SummaryIterable, args...)
    # run iteration for iterable iter
    next = iterate(iter.iter,args...)

    # if iteration happend, display
    if next !== nothing
        iter.summ(next[1])
    end

    # don't do anything to next
    return next
end

summary(iter::I, summ::F) where {I,F} = SummaryIterable{I,F}(iter,summ)


#####
##### History iterable
#####

"""
An iterable to output the history of iterations.
Misfit, iterations taken by the method, and 
stoppign conditions are the ones coming to mind.

A type `IterationHistory` is set as a mutable structure where
a counter of iterations and a `Dict` to keep information generated
during each method's execution. The original implementation, given in

https://github.com/JuliaMath/IterativeSolvers.jl/blob/master/src/history.jl

provides more fields than just

`iters::Int`: number of iterations taken by the method
`data::Dict{Symbol,Any}`: keep info generated by the method.

In my opinion, this is a very smart way of doing things.
A smaller (but not necessarily small) package that handles
a summary of the iterations. We usually check the misfit,
and SNR per iteration, so you should be expecting two Arrays{Float,1}(iters). 
"""

#= 
The iterable, the constructor, the original:
 The original implementation provides other options, such as
 partial summaries. I will just put everything at once, with one
 whole log that will be always part of the output. By this I mean that
 every iterative method being coded in this non-official
 not-registred to-learn-a-lot package should output the solution array
 and a history structure from which we can access number of iterations
 and other measurements of interest.
=#


#=
 For this type, since we will be storing information 
 it will be necessary to define functions like
 setindex and getindex for the custom type.
 This is what you are doing below
=#

# The iterable
mutable struct IterationHistory
    iters::Int
    data::Dict{Symbol, Any}
end

function IterationHistory()
    IterationHistory( 0, Dict{Symbol, Any}() )
end

"""
getindex overload for IterationHistory.

Get collection associated with a given symbol `s` in `IterationHistory`.
"""
getindex(hist::IterationHistory, s::Symbol) = hist.data[s];
getindex(hist::IterationHistory, s::Symbol, kwargs) = hist.data[s][kwargs...];

"""
setindex!(hist, val, s)

setindex overload for IterationHistory.

Set (collection of) elements associated with a given symbol `s` in
`IterationHistory` to val.
"""
setindex!(hist::IterationHistory, val, s::Symbol) = hist.data[s] = val;
setindex!(hist::IterationHistory, val, s::Symbol, kwargs) =  hist.data[s][kwargs...] = val;

"""
push!(hist, key, vec)

Push the contents of `vec` to `key` in `hist`.
`vec` is a collection, here constrained to be either a vector or a
tuple, that contains useful information on your iterations.
"""
function push!(hist::IterationHistory, key::Symbol, vec::Union{Vector, Tuple})
    getkey = hist.data[key] # extract key from dictionary
    shift = size(getkey,2)  # extract its size along dim=2 to shift
    # position in output
    iter = hist.iters       # current iteration
    pos = (iter-1)*shift    # relative position

    for i in 1:min(shift,length(vec))
        getkey[pos+i] = vec[i] # notice how this place vec[i] in hist.data[key][pos+i]
    end
end

push!(hist::IterationHistory, key::Symbol, data) = push_custom_data!(hist, key, data)
push_custom_data!(hist::IterationHistory, key::Symbol, data) = hist.data[key][hist.iters] = data

"""
    reserve!(hist,key,maxIter)
    reserve!(type,hist,key,maxIter)

    Function to reserve space in hist to save the summary.
    It is, of course, based on the maximum number of iterations, but
    it also requires type definitions for what you want to solve. As
    usual, this probably improves performance. Notice, however, that
    the method might converge before maxIter is reached, in which case
    a solution to shrink the unused space was also presented in `shrink!`.
"""
function reserve!(hist::IterationHistory,keys::Vector{Symbol},kwargs...)
    for key in keys
        reserve!(hist, key, kwargs...)
    end
end

function reserve!(hist::IterationHistory, key::Symbol, kwargs...)
    # Default for not type annotated calls
    reserve!(Float64, hist, key, kwargs...)
end

function reserve!(T::Type, hist::IterationHistory, key::Symbol, kwargs...)
    _reserve!(T, hist, key, kwargs...) # an internal method for type
end                                    # annotated calls

function _reserve!(T::Type, hist::IterationHistory, key::Symbol, dim::Int, kwargs...)
    hist.data[key] = Vector{T}(undef,dim)
end

# For future me:
# We can store intermediate solutions and gradients with smth like
# this 
function _reserve!(T::Type, hist::IterationHistory, key::Symbol,dim1::Int, dim2::Int, kwargs...)
    hist.data[key] = Matrix{T}(undef,dim1,dim2)
end

"""
    shrink!(hist::IterationHistory)

    Shirnks the unused space reserved (by reserve!) to log hist.
"""
function shrink!(hist::IterationHistory)
    for key in keys(hist)
        # extract keys from hist
        # /NB/ Keys are actually only defined in the `Dict` hist.data
        # We'll override the function keys
        getKey = hist.data[key]
        if isa(getKey,Vector)
            resize!(getKey, hist.iters)
        elseif isa(getKey, Matrix)
            hist.data[key] = getKey[1:hist.iters,:]
        end
    end
end

#####
##### Set useful functions 
#####

# Can be used to add more info later
function nextiter!(hist::IterationHistory)
    hist.iters += 1;
end

# To extract keys directly from IterationHistory objects
keys(hist::IterationHistory) = keys(hist.data)

# Prediction quality
function prediction_quality(pred::AbstractArray,ideal::AbstractArray)

    α = dot(ideal,pred) / dot(pred,pred);
    Q = 20 * log( norm(ideal,2) / norm(ideal - α .* pred,2)  )

    return Q
end