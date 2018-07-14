module Expectations

# Load external dependencies. 
using Distributions, QuantEcon, QuadGK

# Load internal files. 
include("types.jl"),
include("lib/discrete.jl"),
include("lib/continuous.jl")


end # module
