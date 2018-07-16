# Types for quadrature algorithms. 
abstract type QuadratureAlgorithm end 
abstract type ExplicitQuadratureAlgorithm end 

# Concrete types for quadrature algorithms. 
struct Gaussian <: QuadratureAlgorithm end # Distribution-family specific quadrature.

struct FiniteDiscrete <: ExplicitQuadratureAlgorithm end # Dot-product basically. 
struct Trapezoidal <: ExplicitQuadratureAlgorithm end # For iterable expectations. 

# Abstract types for expectations. 
abstract type Expectation end # Supports E(f)

# Concrete types for expectations. 

#= For an example of using abstract types named in this way, see: https://github.com/JuliaStats/Distributions.jl/blob/2d98eb6f31e9a92cce416e7391a84cff9bba7292/src/truncate.jl#L1-L10. We define a family of Truncated{blahblahblah} types parametrically, but use the abstract Truncated as a supertype for all Truncated distributions. 
=#

struct IterableExpectation{NT, WT} <: Expectation # Supports E(f), nodes, weights, * 
    nodes::NT 
    weights::WT 
end 