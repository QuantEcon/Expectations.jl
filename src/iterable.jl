#= 
    All iterable expectations.
=#

# Callable behavior for the object. Parameters because we cannot add methods to an abstract type. 
function (e::IterableExpectation{NT, WT})(f::Function; kwargs...) where {NT, WT}
    applicable(f, rand([e.nodes...])) || throw(ArgumentError("The function doesn't accept elements from the distribution's support."))
    return dot(f.(e.nodes), e.weights; kwargs...) 
end 

# Getters for the object. 
nodes(e::IterableExpectation) = e.nodes 
weights(e::IterableExpectation) = e.weights 

# Linear operator behavior (TBD).

#= 
    Discrete iterable expectations. 
=# 

# Constructors for the object. 
function expectation(dist::D, alg::Type{FiniteDiscrete} = FiniteDiscrete; kwargs...) where {D <: DiscreteUnivariateDistribution}
    return _expectation(dist, alg; kwargs...)
end

function _expectation(dist::D, alg::Type{FiniteDiscrete}; kwargs...) where {D <: DiscreteUnivariateDistribution}
    hasfinitesupport(dist) || throw(MethodError("Countably infinite distributions are not supported."))
    ourSupport = support(dist)
    ourWeights = pdf.(dist, support(dist))
    sum(ourWeights) ≈ 1.0 || warn("The distribution supplied is not approximately equal to 1 in mass.")
    return IterableExpectation(ourSupport, ourWeights);
end 

#= 
    Continuous iterable expectations (no nodes supplied.)
=#

# General catchall behavior --> Gauss-Legendre quadrature. 
function expectation(dist::D, alg::Type{<:QuadratureAlgorithm} = Gaussian; kwargs...) where {D <: ContinuousUnivariateDistribution}
    return _expectation(dist, alg; kwargs...)
end 

function _expectation(dist::D, alg::Type{Gaussian}; n = 400, kwargs...) where {D <: ContinuousUnivariateDistribution}
    a = minimum(dist)
    b = maximum(dist)
    (a > -Inf && b < Inf) || throw(MethodError("The distribution must be defined on a compact interval."))
    nodes, weights = qnwlege(n, a, b)
    weights = [weights[i] * pdf(dist, nodes[i]) for i in 1:length(nodes)]
    return IterableExpectation(nodes, weights);
end 

# Specific method for normal distributions. 
# Number of points was calibrated by trial.
function _expectation(dist::D, alg::Type{Gaussian}; n = 24, kwargs...) where {D <: Normal}
    μ = mean(dist)
    σ_2 = var(dist)
    nodes, weights = qnwnorm(n, μ, σ_2)
    return IterableExpectation(nodes, weights)
end 

# Specific method for lognormal distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 23, kwargs...) where {D <: LogNormal}
end 

#= 
    Continuous iterable distributions (nodes supplied.)
=#

function expectation(dist::D, nodes::NT, alg::Type{Gaussian}; kwargs...) where {D <: ContinuousUnivariateDistribution, NT}
    return _expectation(dist, nodes, alg; kwargs...)
end

# Put stuff here. 