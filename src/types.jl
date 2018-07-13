"""
A common supertype for all objects in this package. 
"""
abstract type AbstractExpectation end 

"""
Expectation object for 
"""
struct CallableExpectation{T <: Distribution} <: AbstractExpectation
    D::T # The underlying distribution goes here. 
end 

