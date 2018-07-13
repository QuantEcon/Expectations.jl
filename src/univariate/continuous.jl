
# Helper functions 
_expectation(E::AbstractExpectation{<:ContinuousUnivariateDistribution}, alg::Type{Quadrature1}; N=10, kwargs...) = (1,N) #i.e. algorithm 1, N=N
_expectation(E::AbstractExpectation{<:Normal}, alg::Type{Quadrature2}; N=5, kwargs...) = (2,N) #i.e. algorithm 2 specialization, N=N

# expectation() methods
function expectation(E::AbstractExpectation{<:ContinuousUnivariateDistribution}, alg::Type{<:QuadratureAlgorithm} = Quadrature1; kwargs...)
    return _expectation(E, alg; kwargs...)
end 
function expectation(E::AbstractExpectation{<:Normal}, alg::Type{<:QuadratureAlgorithm} = Quadrature2; kwargs...)
    return _expectation(E, alg; kwargs...)
end 