# __precompile__(true)

module Expectations

# Load external dependencies. 
using Compat, FastGaussQuadrature
using Compat.LinearAlgebra
using SpecialFunctions
using Reexport
@reexport using Distributions

# Load internal files. 
include("types.jl"),
include("iterable.jl")

# Export 
export expectation, _expectation, Expectation, IterableExpectation, Gaussian, Trapezoidal, FiniteDiscrete, QuadratureAlgorithm, ExplicitQuadratureAlgorithm, nodes, weights

end # module
