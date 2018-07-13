using Expectations, Distributions

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
v1 = Expectation(Normal(0, 1))
v2 = Expectation(Uniform(0, 1))

@test expectation(v1) == expectation(v1, Quadrature2, N=5) #i.e. defaults to Algorithm 2 with N=5 as the default
@test expectation(v1, N = 7) == expectation(v1, Quadrature2, N=7) #Can change the default value
@test expectation(v1, Quadrature1) == expectation(v1, Quadrature1, N=10) #can use use Quadrature1 (with the different default)
@test expectation(v2) == expectation(v2, Quadrature1, N=10) #i.e. uses algorithm 1 by default
@test expectation(v2, N = 8) == expectation(v2, Quadrature1, N=8) #Can change the default value
@test_throws MethodError  expectation(v2, Quadrature2) #Not defined! Should throw