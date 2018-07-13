# Generic fallback 
function (E::Expectation{T})(f::Function = identity) where {T <: UnivariateDistribution}
    # Could implement some generic LLN type logic here. 
end 

# Univariate distributions 
"""
Method for generic discrete univariate distributions (finite or countable).
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: DiscreteUnivariateDistribution}
    finite = hasfinitesupport(E.D)
    if finite
        distsupport = support(E.D)
        # Check that f is defined on the support. 
        # method_exists(f, eltype(distsupport)) ? true : MethodError("The function is not defined on the distribution's support")
        probs = pdf.(E.D, support(E.D))
        vals = f.(support(E.D))
        # Check that the dot product is defined.  
        # method_exists(dot, (typeof(probs), typeof(vals))) ? true : MethodError("The function values and probabilities are not dottable.")
        return dot(probs, vals)
    else 
        # Use Lazy.jl? Can introduce quadrature rule later. 
        # Should probably call out to a special quadrature function. 
        # Use `eltype` to get even more genericity out of this stuff. 
    end 
end 

"""
Method for generic continuous univariate distributions. 
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: ContinuousUnivariateDistribution}
end 

# Specific univariate distributions. Implements specific quadrature rules from Maranda-Fackler. 
"""
Method for univariate normal distribution. 
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: Normal}
end

"""
Method for univariate lognormal distribution. 
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: LogNormal}
end

"""
Method for univariate (continuous) uniform distribution.
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: Uniform}
    D = E.D
    # Checks
    @assert -Inf < D.a < Inf && -Inf < D.b < Inf "One or more bounds is infinite" # Technically, this isn't necessary, since typeof(Inf) <: Float64. But it's a good check.
    @assert applicable(f, D.a) && applicable(*, 0.5, D.a) && applicable(sum, D.a) "f does not take floats into a space that can be dotted with probabilities"
    # Compute and return.
    return univariate_quad(D, f; kwargs...)
end

"""
Method for univariate (discrete) uniform distribution.
"""
function (E::Expectation{T})(f::Function = identity; kwargs...) where {T <: DiscreteUniform}
    D = E.D   
    # Checks
    @assert -Inf < D.a < Inf && -Inf < D.b < Inf "One or more bounds is infinite"
    @assert applicable(f, D.a) && applicable(*, 0.2, f(D.a)) && applicable(sum, f(D.a)) "f does not take integers into a space that can be dotted with probabilities"
    # Compute and return 
    vals = f.(collect(D.a:1:D.b))
    return D.pv * sum(vals)
end