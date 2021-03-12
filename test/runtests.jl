using Test

# Using Testing Module

println("Testing")
@test true                      #> Test Passed
@test [1, 2] + [2, 1] == [3, 3] #> Test Passed
println("... tests passed!")

include("test_quater.jl")
include("test_qfft.jl")
