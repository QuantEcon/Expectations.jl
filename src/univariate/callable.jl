# Define the behavior of an Expectation as a callable object. 
include("../types.jl")

""" 
Generic function for calling an expectation of a function
"""
function (E::Expectation{T})(f::Function) where {T}
end 

# Univariate Distributions 
"""
Method for generic discrete univariate distributions (finite or countable).
"""
function (E::Expectation{T})(f::Function = x -> x; tol = 1e-8) where {T <: DiscreteUnivariateDistribution}
    # Switch on discreteness.  
    finite = hasfinitesupport(E.D)
    if finite
        distsupport = support(E.D)
        # Check that f is defined on the support. 
        method_exists(f, eltype(distsupport)) ? true : MethodError("The function is not defined on the distribution's support")
        probs = pdf.(E.D, support(E.D))
        vals = f.(support(E.D))
        # Check that the dot product is defined.  
        method_exists(dot, (typeof(probs), typeof(vals))) ? true : MethodError("The function values and probabilities are not dottable.")
        return dot(probs, vals)
    else 
        # Can revise this quadrature rule later. 
        val = 0 
        error = Inf 
        while error > tol 
        end 
    end 
end 