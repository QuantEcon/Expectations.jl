# High-level types. Parametric lets us do stuff like AbstractExpectation{<:ContinuousUnivariateDistribution}.
"""
A common supertype for all expectation-takers in this package. 
"""
abstract type AbstractExpectation{T} end 

"""
A common supertype for all quadrature algorithms in this package. 
"""
abstract type QuadratureAlgorithm end 

"""
Expectation object for callable expectations.
"""
struct NodelessExpectation{T <: Distribution} <: AbstractExpectation{T}
    D::T # The underlying distribution goes here. 
end 

# Type Alias 
const Expectation = NodelessExpectation

# """
# Expectation object for when we fix the nodes. 
# """
# struct NodesExpectation{T <: Distribution, CN <: AbstractArray} where {T <: Distribution, CN <: AbstractArray} <: AbstractExpectation{T} 
#     D::T
#     nodes::CN # AbstractArray supers things like Vector, UnitRange, StepRangeLen, etc. 
    
#     # Generic inner constructor. 
#     NodesExpectation(D, nodes) = eltype(nodes) <: eltype(support(D)) ? new(D, nodes) : ArgumentError("The container does not hold the same type of elements as the distribution support. ")
# end 

# Specific algorithm types. 
abstract type Quadrature1 <: QuadratureAlgorithm end 
abstract type Quadrature2 <: QuadratureAlgorithm end 