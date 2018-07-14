# Generic getters. 
nodes(e::AbstractExpectation) = e.nodes 
weights(e::AbstractExpectation) = e.weights 


# Get this working once https://github.com/JuliaLang/julia/issues/14919 is fixed. 
# # Generic evaluator. 
# function (e::AbstractExpectation)(f::Function = identity) 
#     @assert applicable(f, rand(e.nodes)) "The function f does not accept arguments from the sample space."
#     return e.func(f)
# end 

