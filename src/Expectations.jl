# __precompile__(true)

module Expectations

# Load external dependencies.
using FastGaussQuadrature
using LinearAlgebra
using SpecialFunctions
using Distributions

# Load internal files.
include("types.jl"),
include("iterable.jl")
include("mixturemodels.jl")

# Export
export expectation, Expectation, IterableExpectation, Gaussian, Trapezoidal, FiniteDiscrete, QuadratureAlgorithm, ExplicitQuadratureAlgorithm, nodes, weights, QuantileRange
export MixtureExpectation, expectations

end # module
