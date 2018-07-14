#= 
    Generic non-linear methods (i.e., those without user-supplied nodes.)
=#

# Methods for continuous distributions under AdaptiveGaussian quadrature. 
function _expectation(dist::D, alg::Q; n = 200, kwargs...) where {D <: ContinuousDistribution, Q <: AdaptiveGaussian}
    @assert length(dist) == 1 "The QuadGK implementation doesn't yet support integrating over high-dimensional spaces."
    # Get the nodes and weights.
    nodes, weights = konrod(n) # TEST THIS. Add Gaussian? 
    # Construct the functor. 
    a, b = extrema(support(dist))
    ourPDF = x -> pdf(dist, x)
    func = f -> quadgk(x -> f(x) * ourPDF(x), a, b; kwargs...)
    # Write the algorithm. 
    return Continuous(nodes, weights, func, alg)
end

# MethodError for unsupported distributions under specialized quadrature 
function _expectation(dist::D, alg::Q; kwargs...) where {D <: ContinuousDistribution, Q <: Gaussian}
    throw(ArgumentError("Specialized quadrature is not yet available for this distribution."))
end 

# MethodError for continuous distributions under FiniteDiscrete quadrature. 
function _expectation(dist::D, alg::Q; kwargs...) where {D <: ContinuousDistribution, Q <: FiniteDiscrete}
    throw(ArgumentError("The Finite Difference algorithm is not supported for general continuous distributions."))
end 

# Method to make a generic continuous expectation object given a continuous univariate distribution
function expectation(dist::D, alg::Q = AdaptiveGaussian; kwargs...) where {D <: ContinuousDistribution, Q <: QuadratureAlgorithm}
    return _expectation(dist, alg; kwargs...)
end 

#= 
    Evaluator. Will genericize this once https://github.com/JuliaLang/julia/issues/14919 is working. 
=#

function (e::Continuous)(f::Function = identity) 
    @assert applicable(f, rand(e.nodes)) "The function f does not accept arguments from the sample space."
    return e.func(f)
end 

#= 
    Specialized non-linear methods.
=#

# Normal distribution, using the QuantEcon method. 
function _expectation(dist::D, alg::Q; n = 200, kwargs...) where {D <: Normal, Q <: Gaussian}
    μ = mean(dist)
    σ = std(dist)
    nodes, weights = qnwnorm(n, μ, σ^2)
    func = f -> dot(f.(nodes), weights)
    return Continuous(nodes, weights, func, alg)
end 

# Lognormal distribution, using the QuantEcon method. 
function _expectation(dist::D, alg::Q; n = 200, kwargs...) where {D <: LogNormal, Q <: Gaussian}
    μ = mean(dist)
    σ = std(dist)
    nodes, weights = qnwlogn(n, μ, σ^2)
    func = f -> dot(f.(nodes), weights)
    return Continuous(nodes, weights, func, alg)
end 

# etc.

#= 
    Generic linear methods (i.e., those with nodes supplied.)
=#

# Methods for general linear expectations under Trapezoidal quadrature. 
function _expectation(dist::D, nodes, alg::Q; n = 200, kwargs...) where {D <: ContinuousDistribution, Q <: Trapezoid}
    @assert length(dist) == 1 "For now, only univariate distributions are supported."

    # Compute weights. 
    M = length(nodes)
    Δ = diff(nodes)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1.
    ourPDF = x -> pdf(dist, x)
    f_vec = ourPDF.(nodes)
    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    weights = [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]

    # Compute other objects and return. 
    ourFunc = f -> dot(f.(nodes), weights)
    @assert ourFunc(identity) ≈ 1 atol = 1e-2 "The grid doesn't capture enough probability mass."
    return LinearContinuous(nodes, weights, ourFunc, alg)
end 

# Method to make a generic linear continuous expectation object given a continuous univariate distribution and nodes
function expectation(dist::D, nodes::AbstractArray, alg::Q = Trapezoid; kwargs...) where {D <: ContinuousDistribution, Q <: QuadratureAlgorithm}
    return _expectation(dist, nodes, alg; kwargs...)
end 
