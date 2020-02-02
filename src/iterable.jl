#=
    All iterable expectations.
=#

# Callable behavior for the object. Parameters because we cannot add methods to an abstract type.
"""
    function (e::IterableExpectation{NT, WT})(f::Function; kwargs...) where {NT, WT}

Implements callable behavior for `IterableExpectation` objects.
"""
function (e::IterableExpectation{NT, WT})(f::Function; kwargs...) where {NT, WT}
    applicable(f, rand(nodes(e))) || throw(ArgumentError("The function doesn't accept elements from the distribution's support."))
    return dot(f.(nodes(e)), weights(e))
end

# Getters for the object.
"""
    nodes(e::IterableExpectation)

Returns the nodes of an `IterableExpectation`.
"""
nodes(e::IterableExpectation) = e.nodes

"""
    weights(e::IterableExpectation)

Returns the weights of an `IterableExpectation`.
"""
weights(e::IterableExpectation) = e.weights

# Linear operator behavior.
import Base.*

# Right-multiplying an expectation by something.
"""
    *(e::IterableExpectation, h::AbstractArray) = dot(h, weights(e))

Implements the right-application of an `IterableExpectation` by a vector of values on its nodes.
"""
*(e::IterableExpectation, h::AbstractArray) = dot(h, weights(e))


# Left-multiplying an expectation by a scalar.
"""
    *(r::Real, e::IterableExpectation) =  IterableExpectation(nodes(e), r * weights(e))

Implements left-multiplication of an `IterableExpectation` by a real scalar.
"""
*(r::Real, e::IterableExpectation) =  IterableExpectation(nodes(e), r * weights(e)) # Necessary because, for example, multiplying UnitRange * 2 = StepRange


#=
    Discrete iterable expectations.
=#

# Constructors for the object.
"""
    expectation(dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete} = FiniteDiscrete; kwargs...) = _expectation(dist, alg; kwargs...)

Dispatcher for (finite) discrete univariate expectations.
"""
expectation(dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete} = FiniteDiscrete; kwargs...) = _expectation(dist, alg; kwargs...)

"""
    function _expectation(dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete}; kwargs...)

Auxiliary constructor for an `IterableExpectation` object.
"""
function _expectation(dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete}; kwargs...)
    hasfinitesupport(dist) || throw(ArgumentError("Countably infinite distributions are not supported."))
    ourSupport = support(dist)
    ourWeights = pdf.(Ref(dist), support(dist))
    sum(ourWeights) ≈ 1.0 || warn("The distribution supplied is not approximately equal to 1 in mass.")
    return IterableExpectation(ourSupport, ourWeights);
end

#=
    Continuous iterable expectations (no nodes supplied.)
=#

# General catchall behavior --> Gauss-Legendre quadrature.
"""
    expectation(dist::ContinuousUnivariateDistribution, alg::Type{<:QuadratureAlgorithm} = Gaussian; kwargs...) = _expectation(dist, alg; kwargs...)

Dispatcher for continuous univariate expectations.
"""
expectation(dist::ContinuousUnivariateDistribution, alg::Type{<:QuadratureAlgorithm} = Gaussian; kwargs...) = _expectation(dist, alg; kwargs...)

"""
    function _expectation(dist::ContinuousUnivariateDistribution, alg::Type{Gaussian}; n = 500, kwargs...)

Implements Gauss-Legendre quadrature for continuous univariate distributions for which no specialized method exists.
"""
function _expectation(dist::ContinuousUnivariateDistribution, alg::Type{Gaussian}; n = 500, kwargs...)
    a = minimum(dist)
    b = maximum(dist)
    (a > -Inf && b < Inf) || throw(ArgumentError("The distribution must be defined on a compact interval."))
    rawNodes, rawWeights = gausslegendre(n)
    # Transform nodes to proper interval.
    nodes = map(x -> (0.5(b-a))*x + (a+b)/2, rawNodes)
    # Add pdf to weights.
    compoundWeights = [rawWeights[i] * pdf(dist, nodes[i]) for i in 1:n]
    # Add scale factor to weights.
    weights = (b-a)/2 * compoundWeights
    return IterableExpectation(nodes, weights);
end

# QNWDist implementation.
"""
    function _expectation(dist::ContinuousUnivariateDistribution, alg::Type{QuantileRange}; n::Int = 50, q0::Real = 0.001, qN::Real = 0.999 kwargs...)

Implementation of the qnwdist() quadrature scheme written by Spencer Lyon (PhD. NYU), as part of the QuantEcon.jl library. Used with permission. For further details, see: https://github.com/QuantEcon/QuantEcon.jl/blob/be0a32ec17d1f5b04ed8f2e52604c70c69f416b2/src/quad.jl#L892.
"""
function _expectation(dist::ContinuousUnivariateDistribution, alg::Type{QuantileRange}; n::Int = 50, q0::Real = 0.001, qN::Real = 0.999, kwargs...)
    # check nondegeneracy.
    all(isfinite.(params(dist))) || throw(ArgumentError("Distribution must be nondegenerate."))
    # _quadnodes in the QuantEcon.
    quantiles = Compat.range(q0, stop=qN, length=n)
    nodes = quantile.(Ref(dist), quantiles)
    # qnwdist in the QuantEcon.
    weights = zeros(n)
    for i in 2:n-1
        weights[i] = cdf(dist, (nodes[i] + nodes[i+1])/2) - cdf(dist, (nodes[i] + nodes[i-1])/2)
    end
    weights[1] = cdf(dist, (nodes[1] + nodes[2])/2)
    weights[end] = 1 - cdf(dist, (nodes[end-1] + nodes[end])/2)
    return IterableExpectation(nodes, weights)
end

# Specific method for normal distributions.
# Number of points was calibrated by trial.
"""
    function _expectation(dist::Normal, alg::Type{Gaussian}; n = 30, kwargs...)

Implements Gauss-Hermite quadrature for normal distributions.
"""
function _expectation(dist::Normal, alg::Type{Gaussian}; n = 30, kwargs...)
    σ = std(dist)
    μ = mean(dist)
    (isfinite(σ) && isfinite(μ)) || throw(ArgumentError("Parameters σ, μ must be finite."))
    gh = gausshermite(n)
    nodes = gh[1].*(sqrt(2)*(σ)) .+ μ
    weights = gh[2]./sqrt(pi)
    return IterableExpectation(nodes, weights)
end

# Specific method for lognormal distributions.
"""
    function _expectation(dist::LogNormal, alg::Type{Gaussian}; n = 30, kwargs...)

Implements Gauss-Hermite quadrature for lognormal distributions.
"""
function _expectation(dist::LogNormal, alg::Type{Gaussian}; n = 30, kwargs...) # Same settings for the normal method.
    m = mean(dist)
    v = var(dist)
    (isfinite(m) && isfinite(v)) || throw(ArgumentError("Infinite μ or σ^2 are not supported."))
    # get normal nodes
    gh = gausshermite(n)
    μ = log(m^2/sqrt(v + m^2))
    σ = sqrt(log(v/m^2 + 1))
    nodes = gh[1].*sqrt(2).*(σ) .+ μ
    weights = gh[2]./sqrt(pi)
    # get new nodes
    map!(x -> exp(x), nodes, nodes)
    return IterableExpectation(nodes, weights) # Transform the output.
end

# Specific method for beta distributions.
"""
    function _expectation(dist::Beta, alg::Type{Gaussian}; n = 32, kwargs...)

Implements Gauss-Jacobi quadrature for beta distributions.
"""
function _expectation(dist::Beta, alg::Type{Gaussian}; n = 32, kwargs...)
    α, β = params(dist)
    (isfinite(α) && isfinite(β)) || throw(ArgumentError("The beta distribution supplied is malformed."))
    gj = FastGaussQuadrature.JacobiRec(n, α-1, β-1)
    G = gamma(α)*gamma(β)/gamma(α+β)
    nodes = (1 .- gj[1])/2
    weights = gj[2]/((2.0^(α+β-1.0))*G)
    return IterableExpectation(nodes, weights)
end

# Specific method for exponential distributions.
"""
    function _expectation(dist::Exponential, alg::Type{Gaussian}; n = 32, kwargs...)

Implements Gauss-Laguerre quadrature for Exponential distributions.
"""
function _expectation(dist::Exponential, alg::Type{Gaussian}; n = 32, kwargs...)
    θ = inv(dist.θ) # To correct for the Distributions parametrization.
    isfinite(θ) || throw(ArgumentError("The beta distribution supplied is malformed."))
    gl = gausslaguerre(n)
    nodes = gl[1]./θ
    weights = gl[2]
    return IterableExpectation(nodes, weights)
end

# Specific method for gamma distributions.
"""
    function _expectation(dist::Union{Gamma,Erlang}, alg::Type{Gaussian}; n = 32, kwargs...)

Implements Gauss-Laguerre quadrature for Gamma distributions.
"""
function _expectation(dist::Union{Gamma,Erlang}, alg::Type{Gaussian}; n = 32, kwargs...)
    α, θ = params(dist)
    (isfinite(θ) && isfinite(θ)) || throw(ArgumentError("The beta distribution supplied is malformed."))
    gl = gausslaguerre(n, α-1)
    nodes = gl[1]./θ
    weights = gl[2]./gamma(α)
    return IterableExpectation(nodes, weights)
end

#=
    Continuous iterable distributions (nodes supplied.)
=#

# Dispatcher.
"""
    expectation(dist::ContinuousUnivariateDistribution, nodes, alg::Type{<:ExplicitQuadratureAlgorithm} = Trapezoidal; kwargs...)  = _expectation(dist, nodes, alg; kwargs...)

Dispatcher for distributions with user-defined nodes.
"""
expectation(dist::ContinuousUnivariateDistribution, nodes, alg::Type{<:ExplicitQuadratureAlgorithm} = Trapezoidal; kwargs...)  = _expectation(dist, nodes, alg; kwargs...)

# Trapezoidal general behavior.
"""
    function _expectation(dist, nodes::AbstractArray, alg::Type{Trapezoidal}; kwargs...)

Implements trapezoidal integration for general distributions with user-defined nodes.
"""
function _expectation(dist, nodes::AbstractArray, alg::Type{Trapezoidal}; kwargs...)
    isfinite(minimum(dist)) && isfinite(maximum(dist)) || throw(ArgumentError("The distribution must be defined on a compact interval."))
    (first(nodes) >= minimum(dist)) && (last(nodes) <= maximum(dist)) || throw(ArgumentError("The nodes exceed the distribution's support."))
    M = length(nodes)
    Δ = diff(nodes)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1. And NaN will throw an error if we try to use it.
    f_vec = pdf.(Ref(dist), nodes)
    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    allWeights = [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]
    return IterableExpectation(nodes, allWeights)
end

# Trapezoidal for regular.
"""
    @compat function _expectation(dist, nodes::AbstractRange, alg::Type{Trapezoidal}; kwargs...)

Overloads trapezoidal integration for cases where the user-defined grids are regular.
"""
@compat function _expectation(dist, nodes::AbstractRange, alg::Type{Trapezoidal}; kwargs...)
    isfinite(minimum(dist)) && isfinite(maximum(dist)) || throw(ArgumentError("The distribution must be defined on a compact interval."))
    (first(nodes) >= minimum(dist)) && (last(nodes) <= maximum(dist)) || throw(ArgumentError("The nodes exceed the distribution's support."))
    M = length(nodes)
    Δ = nodes[2] - nodes[1]
    f_vec = pdf.(Ref(dist), nodes)
    interiorWeights = [f_vec[i] * Δ for i = 2:M-1]
    allWeights = [f_vec[1]/2 * Δ; interiorWeights; f_vec[M]/2 * Δ]
    return IterableExpectation(nodes, allWeights)
end

#=
    Convenience functions.
=#
"""
    expectation(f::Function, dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete} = FiniteDiscrete; kwargs...) = expectation(dist, alg; kwargs...)(f)

Convenience function for (finite) discrete univariate distributions.
"""
expectation(f::Function, dist::DiscreteUnivariateDistribution, alg::Type{FiniteDiscrete} = FiniteDiscrete; kwargs...) = expectation(dist, alg; kwargs...)(f)

"""
    expectation(f::Function, dist::ContinuousUnivariateDistribution, alg::Type{<:QuadratureAlgorithm} = Gaussian; kwargs...) = expectation(dist, alg; kwargs...)(f)

Convenience function for continuous univariate distributions.
"""
expectation(f::Function, dist::ContinuousUnivariateDistribution, alg::Type{<:QuadratureAlgorithm} = Gaussian; kwargs...) = expectation(dist, alg; kwargs...)(f)

"""
    expectation(f::Function, dist::ContinuousUnivariateDistribution, nodes::AbstractArray, alg::Type{<:ExplicitQuadratureAlgorithm} = Trapezoidal; kwargs...) = expectation(dist, nodes, alg; kwargs...)(f)

Convenience function for continuous univariate distributions with user-supplied nodes.
"""
expectation(f::Function, dist::ContinuousUnivariateDistribution, nodes::AbstractArray, alg::Type{<:ExplicitQuadratureAlgorithm} = Trapezoidal; kwargs...) = expectation(dist, nodes, alg; kwargs...)(f)
