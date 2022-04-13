__precompile__()

module HCDSP

using Statistics, StatsBase, NaNMath
using Distributed, DelimitedFiles
using LinearAlgebra, FFTW
using TSVD, DSP
using Random

# Include quternions
include("./Quaternions/Fundamentals/Quaternions.jl")
include("./Quaternions/QFT.jl")
include("./Quaternions/QSVD.jl")

# include SNR
include("./noise/noise.jl")
include("./denoise/denoise.jl")
include("./operators/operators.jl")
include("./stdlib/stdlib.jl")

# include IterativeMethods
include("./optimization/iterables.jl")
include("./optimization/backtrack.jl")
include("./optimization/solvers/solvers.jl")

greet() = print("Hello World! This is a HyperComplex Digital Signal Processing toolbox for us to do signal processing with quaternion numbers.")

end # module
