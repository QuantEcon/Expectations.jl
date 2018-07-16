# Types for quadrature algorithms. 
abstract type QuadratureAlgorithm end 
abstract type IterableQuadratureAlgorithm end 

# Concrete types for quadrature algorithms. 
struct AdaptiveGaussian <: QuadratureAlgorithm end # Catchall. 
struct Gaussian <: QuadratureAlgorithm end # Distribution-family specific quadrature.
struct FiniteDiscrete <: QuadratureAlgorithm end # Dot-product basically. 

struct Trapezoidal <: IterableQuadratureAlgorithm end # For Iterable expectations. 

# Abstract types for expectations. 
abstract type AbstractContinuousExpectation end 
abstract type AbstractIterableContinuousExpectation <: AbstractContinuousExpectation end

abstract type AbstractDiscreteExpectation end 

# Concrete types for expectations. 
struct ContinuousUnivariateExpectation <: AbstractContinuousExpectation
    func::Function
end 

struct IterableContinuousUnivariateExpectation{NT, WT} <: AbstractIterableContinuousExpectation
    nodes::NT 
    weights::WT 
end 

struct DiscreteUnivariateExpectation{ST, WT} <: AbstractDiscreteExpectation
    support::ST 
    weights::WT
end 