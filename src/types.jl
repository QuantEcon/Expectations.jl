# Types for quadrature algorithms. 
abstract type QuadratureAlgorithm end 
struct AdaptiveGaussian <: QuadratureAlgorithm end 
struct Gaussian <: QuadratureAlgorithm end
struct FiniteDiscrete <: QuadratureAlgorithm end
struct Trapezoid <: QuadratureAlgorithm end 

# Abstract types for expectations. 
abstract type AbstractExpectation end 

abstract type AbstractContinuous <: AbstractExpectation end # Don't need discrete, since only doing finite. 

# Concrete types for distributions 
struct LinearContinuous <: AbstractContinuous
    nodes::AbstractArray
    weights::AbstractArray
    func::Function
    alg::QuadratureAlgorithm
end 

struct Continuous <: AbstractContinuous
    nodes::AbstractArray 
    weights::AbstractArray
    func::Function
    alg::QuadratureAlgorithm
end


struct Discrete <: AbstractExpectation
    nodes::AbstractArray 
    weights::AbstractArray # Lets us store either a vector of values or a function. 
    func::Function
    alg::QuadratureAlgorithm
end 