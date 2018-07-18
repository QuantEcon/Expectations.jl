#= 
    All iterable expectations.
=#

# Callable behavior for the object. Parameters because we cannot add methods to an abstract type. 
function (e::IterableExpectation{NT, WT})(f::Function; kwargs...) where {NT, WT}
    applicable(f, rand(nodes(e))) || throw(MethodError("The function doesn't accept elements from the distribution's support."))
    return dot(f.(nodes(e)), weights(e))
end 

# Getters for the object. 
nodes(e::IterableExpectation) = e.nodes 
weights(e::IterableExpectation) = e.weights 

# Linear operator behavior.
import Base.*

# Right-multiplying an expectation by something. 
function *(e::IterableExpectation, h::AbstractArray)
    return dot(h, weights(e))
end 

# Left-multiplying an expectation by a scalar.
function *(r::Real, e::IterableExpectation)
    return IterableExpectation(nodes(e), r * weights(e)) # Necessary because, for example, multiplying UnitRange * 2 = StepRange
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
    ourWeights = pdf.(Ref(dist), support(dist))
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
    rawNodes, rawWeights = gausslegendre(n)
    # Transform nodes to proper interval. 
    nodes = map(x -> (0.5(b-a))*x + (a+b)/2, rawNodes)
    # Add pdf to weights. 
    compoundWeights = [rawWeights[i] * pdf(dist, nodes[i]) for i in 1:n]
    # Add scale factor to weights.
    weights = (b-a)/2 * compoundWeights
    return IterableExpectation(nodes, weights);
end

# Specific method for normal distributions. 
# Number of points was calibrated by trial.
function _expectation(dist::D, alg::Type{Gaussian}; n = 30, kwargs...) where {D <: Normal}
    σ = std(dist)
    μ = mean(dist)
    (isfinite(σ) && isfinite(μ)) || throw(MethodError("Parameters σ, μ must be finite."))
    gh = gausshermite(n)
    nodes = gh[1].*(sqrt(2)*(σ + μ))
    weights = gh[2]./sqrt(pi)
    return IterableExpectation(nodes, weights)
end 

# Specific method for lognormal distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 10, kwargs...) where {D <: LogNormal} # Same settings for the normal method.
    m = mean(dist)
    v = var(dist)
    (isfinite(m) && isfinite(v)) || throw(MethodError("Infinite μ or σ^2 are not supported."))
    # get normal nodes
    gh = gausshermite(n)
    μ = log(m^2/sqrt(v + m^2))
    σ = sqrt(log(v/m^2 + 1))
    nodes = gh[1].*(sqrt(2)*(σ + μ))
    weights = gh[2]./sqrt(pi)
    # get new nodes 
    map!(x -> exp(x), nodes, nodes)
    return IterableExpectation(nodes, weights) # Transform the output. 
end 

# Specific method for beta distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 32, kwargs...) where {D <: Beta} 
    α, β = params(dist)
    (isfinite(α) && isfinite(β)) || throw(MethodError("The beta distribution supplied is malformed."))
    gj = FastGaussQuadrature.JacobiRec(n, α-1, β-1)
    G = gamma(α)*gamma(β)/gamma(α+β)
    nodes = (1 .- gj[1])/2
    weights = gj[2]/((2.0^(α+β-1.0))*G)
    return IterableExpectation(nodes, weights)
end 

# Specific method for exponential distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 32, kwargs...) where {D <: Exponential} 
    θ = params(dist)[1]
    isfinite(θ) || throw(MethodError("The beta distribution supplied is malformed."))
    gl = gausslaguerre(n)
    nodes = gl[1]./θ
    weights = gl[2]
    return IterableExpectation(nodes, weights)
end 

# Specific method for gamma distributions. 
function _expectation(dist::D, alg::Type{Gaussian}; n = 32, kwargs...) where {D <: Gamma} 
    α, θ = params(dist)
    (isfinite(θ) && isfinite(θ)) || throw(MethodError("The beta distribution supplied is malformed."))
    gl = gausslaguerre(n, α-1)    
    nodes = gl[1]./θ
    weights = gl[2]./gamma(α)
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
function _expectation(dist, nodes::AbstractArray, alg::Type{Trapezoidal}; kwargs...)
    M = length(nodes)
    Δ = diff(nodes)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1. And NaN will throw an error if we try to use it.
    f_vec = pdf.(Ref(dist), nodes)
    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    allWeights = [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]
    return IterableExpectation(nodes, allWeights)
end 

# Trapezoidal for regular. 
@compat function _expectation(dist, nodes::StepRangeLen, alg::Type{Trapezoidal}; kwargs...)
    (first(nodes) >= minimum(dist) && last(nodes) <= maximum(dist)) || throw(ArgumentError("The nodes exceed the distribution's support."))
    M = length(nodes)
    Δ = nodes[2] - nodes[1]
    f_vec = pdf.(Ref(dist), nodes)
    weights = Δ/2 * [f_vec[i] * ((i > 1 && i < M) ? 2 : 1) for i in 1:M]
    return IterableExpectation(nodes, weights)
end 