module Expectations

# Load external dependencies. 
using Distributions, QuantEcon, QuadGK

# Load internal files. 
include("types.jl")
include("univariate/continuous.jl")

# Export whatever we're exporting
export AbstractExpectation, QuadratureAlgorithm, Expectation, Quadrature1, Quadrature2, expectation
"""
Docstring goes here
"""
Expectations

end # module
