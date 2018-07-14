# Method to create an expectation object from a finite discrete distribution. 
function expectation(dist::D; kwargs...) where {D <: DiscreteDistribution}
    @assert hasfinitesupport(dist) "Countably infinite distributions are not currently supported."
    ourPDF = x -> pdf(dist, x)
    ourSupport = support(dist)
    ourFunctional = f -> dot(ourPDF.(ourSupport), f.(ourSupport))
    return Discrete(ourSupport, ourPDF.(ourSupport), ourFunctional, FiniteDiscrete)
end 

# Evaluator. Will genericize this once https://github.com/JuliaLang/julia/issues/14919 is working. 
# Generic evaluator. 
function (e::Discrete)(f::Function = identity) 
    @assert applicable(f, rand(e.nodes)) "The function f does not accept arguments from the sample space."
    return e.func(f)
end 

