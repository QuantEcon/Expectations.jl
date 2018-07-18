# __precompile__(true)

module Expectations

# Load external dependencies. 
using Distributions, Compat, FastGaussQuadrature
using Compat.LinearAlgebra
using Compat.MathConstants

# Load internal files. 
include("types.jl"),
include("iterable.jl")

# Export 
export expectation, _expectation, Expectation, IterableExpectation, Gaussian, Trapezoidal, FiniteDiscrete, QuadratureAlgorithm, ExplicitQuadratureAlgorithm, nodes, weights

end # module
