module Expectations

# Load external dependencies. 
using Distributions, QuantEcon, QuadGK

# Load internal files. 
include("types.jl"),
include("iterable.jl")

# Export 
export expectation, Expectation, IterableExpectation

end # module
