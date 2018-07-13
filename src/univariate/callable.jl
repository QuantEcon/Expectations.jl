# Define the behavior of an Expectation as a callable object. 
include("../types.jl")

# Generic fallback 
function (E::Expectation{T})(f::Function = x -> x) where {T <: UnivariateDistribution} 
    # Could implement some generic LLN type logic here. 
end 

# Univariate distributions 
"""
Method for generic discrete univariate distributions (finite or countable).
"""
function (E::Expectation{T})(f::Function = x -> x; tol = 1e-8) where {T <: DiscreteUnivariateDistribution}
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
    end 
end 

"""
Method for generic continuous univariate distributions.
"""
function (E::Expectation{T})(f::Function = x -> x; tol = 1e-8) where {T <: ContinuousUnivariateDistribution}
end 

# Specific univariate distributions