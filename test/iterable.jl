#= 
    Tests for some discrete distributions
=# 

# Test callable behavior for discrete uniform distribution (support 1:10) 
testDist = DiscreteUniform(1, 10)
E_1 = expectation(testDist)
@test E_1(x -> x) ≈ 5.5
@test E_1(x -> x^2) ≈ 38.5

# Test for error handling. 
testDist2 = Poisson(3)
@test_throws MethodError expectation(testDist2) # Catches unbounded distributions. 
@test_throws MethodError E_1(x -> dot(x, ones(7))) # Test for non-applicable functions. 

#= 
    Tests for some continuous distributions (no nodes supplied.)
=#

# Test callable behavior for truncated normal. (Gauss-Legendre)
testDist3 = Truncated(Normal(0.5), 0, 10)
E_3 = expectation(testDist3)
@test E_3(x -> x) ≈ 1.0091604338370335

# Test for normal. 
testDist4 = Normal()
E_4 = expectation(testDist4)
@test_broken E_4(x -> abs(x)) ≈ sqrt(2/pi) atol = 1e-15 # Most of the weights come up 0, even for small numbers of points. 

#= 
    Tests for some continuous distributions (nodes supplied.) 
=#