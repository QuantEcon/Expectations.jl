# __precompile__(true)

module Expectations

# Load external dependencies. 
using Compat, FastGaussQuadrature
using Compat.LinearAlgebra
import Compat.range
using SpecialFunctions
using Distributions

# Load internal files. 
include("types.jl"),
include("iterable.jl")

# Export 
export expectation, Expectation, IterableExpectation, Gaussian, Trapezoidal, FiniteDiscrete, QuadratureAlgorithm, ExplicitQuadratureAlgorithm, nodes, weights, QuantileRange

end # module
