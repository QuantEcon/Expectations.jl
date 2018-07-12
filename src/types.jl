"""
A common supertype for all expectation operators.
"""
abstract type AbstractExpectation end 

"""
Parametric Expectation object for Distributions. 
"""
struct Expectation{T <: Distribution} <: AbstractExpectation
    D::T # The underlying distribution goes here. 
end 