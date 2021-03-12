__precompile__()

module HCDSP

using LinearAlgebra, FFTW
using Random

# Include quternions
include("./Fundamentals/Quaternions.jl")
include("QFT.jl")

greet() = print("Hello World! This is a HyperComplex Digital Signal Processing toolbox for us to do signal processing with quaternion numbers.")

end # module
