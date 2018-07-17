#= 
    All iterable expectations.
=#

# Callable behavior for the object. Parameters because we cannot add methods to an abstract type. 
function (e::IterableExpectation{NT, WT})(f::Function; kwargs...) where {NT, WT}
    applicable(f, rand([e.nodes...])) || throw(MethodError("The function doesn't accept elements from the distribution's support."))
    return f.(e.nodes)' * e.weights 
end 

# Getters for the object. 
nodes(e::IterableExpectation) = e.nodes 
weights(e::IterableExpectation) = e.weights 

# Linear operator behavior.
import Base.*

# Right-multiplying an expectation by something. 
function *(e::IterableExpectation, h::AbstractArray)
    return h' * weights(e)
end 

# Left-multiplying an expectation by a scalar.
function *(r::Real, e::IterableExpectation)
    return IterableExpectation(e.nodes, r * e.weights) # Necessary because, for example, multiplying UnitRange * 2 = StepRange
end

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
function _expectation(dist::D, alg::Type{Gaussian}; n = 30, kwargs...) where {D <: Normal}
    μ = mean(dist)
    σ_2 = var(dist)
    (isfinite(μ) && isfinite(σ_2)) || throw(MethodError("Infinite μ or σ^2 are not supported."))
    nodes, weights = qnwnorm(n, μ, σ_2)
    return IterableExpectation(nodes, weights)
end 

# Specific method for lognormal distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 100, kwargs...) where {D <: LogNormal} # Same settings for the normal method, since this calls qnwnorm in QuantEcon. 
    μ = mean(dist)
    σ_2 = var(dist)
    (isfinite(μ) && isfinite(σ_2)) || throw(MethodError("Infinite μ or σ^2 are not supported."))
    nodes, weights = qnwlogn(n, μ, σ_2)
    return IterableExpectation(log.(nodes), weights) # Transform the output. 
end 

# Specific method for beta distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 20, kwargs...) where {D <: Beta} # Same settings for the normal method, since this calls qnwnorm in QuantEcon. 
    α, β = params(dist)
    (isfinite(α) && isfinite(β)) || throw(MethodError("The beta distribution supplied is malformed."))
    nodes, weights = qnwbeta(n, α, β)
    return IterableExpectation(nodes, weights)
end 



#= 
    Continuous iterable distributions (nodes supplied.)
=#

# Dispatcher.
function expectation(dist::D, nodes::NT, alg::Type{<:ExplicitQuadratureAlgorithm} = Trapezoidal; kwargs...) where {D <: ContinuousUnivariateDistribution, NT}
    return _expectation(dist, nodes, alg; kwargs...)
end

# Trapezoidal general behavior. 
function _expectation(dist, nodes, alg::Type{Trapezoidal}; kwargs...)
    M = length(nodes)
    Δ = diff(nodes)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1. And NaN will throw an error if we try to use it.
    f = x -> pdf(dist, x)
    f_vec = f.(nodes)
    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    allWeights = [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]
    return IterableExpectation(nodes, allWeights)
end 
