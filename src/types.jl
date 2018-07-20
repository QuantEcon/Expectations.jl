# Types for quadrature algorithms.
"""
Abstract type for quadrature algorithms without user-defined nodes (e.g., Gaussian quadrature.)
"""
abstract type QuadratureAlgorithm end 

"""
Abstract type for quadrature algorithms with user-defined nodes (e.g., trapezoidal integration).
"""
abstract type ExplicitQuadratureAlgorithm end 

# Concrete types for quadrature algorithms.
"""
Gaussian quadrature. See specific methods for what precise algorithm is used (e.g., Gauss-Legendre, Gauss-Hermite, etc.) 
"""
struct Gaussian <: QuadratureAlgorithm end # Distribution-family specific quadrature.

"""
A dot product of a (finite) PDF vector and a finite set of transformed nodes. 
"""
struct FiniteDiscrete <: ExplicitQuadratureAlgorithm end # Dot-product basically.

"""
Trapezoidal integration.
"""
struct Trapezoidal <: ExplicitQuadratureAlgorithm end # For iterable expectations. 

# Abstract types for expectations. 
"""
Abstract type for all expectations.
"""
abstract type Expectation end # Supports E(f)

# Concrete types for expectations. 

#= For an example of using abstract types named in this way, see: https://github.com/JuliaStats/Distributions.jl/blob/2d98eb6f31e9a92cce416e7391a84cff9bba7292/src/truncate.jl#L1-L10. We define a family of Truncated{blahblahblah} types parametrically, but use the abstract Truncated as a supertype for all Truncated distributions. 
=#
"""
Expectations which are paramterized by a vector of nodes (e.g., a discretized support) and corresponding quadrature weights.
"""
struct IterableExpectation{NT, WT} <: Expectation # Supports E(f), nodes, weights, * 
    nodes::NT 
    weights::WT 
end 