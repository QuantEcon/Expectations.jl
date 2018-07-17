# __precompile__(true)

module Expectations

# Load external dependencies. 
using Distributions, QuantEcon, QuadGK

# Load internal files. 
include("types.jl"),
include("iterable.jl")

# Export 
export expectation, _expectation, Expectation, IterableExpectation, Gaussian, Trapezoidal, FiniteDiscrete, QuadratureAlgorithm, ExplicitQuadratureAlgorithm, nodes, weights

end # module
