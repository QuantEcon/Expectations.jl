#= 
    Tests for some discrete distributions
=# 

# Test callable behavior for discrete uniform distribution (support 1:10) 
testDist = DiscreteUniform(1, 10)
E_1 = expectation(testDist)
@test E_1(x -> x) ≈ 5.5
@test E_1(x -> x^2) ≈ 38.5

# Test for linear operator behavior. 
h(x) = 2*x
z = nodes(E_1)
@test E_1 * h.(z) == E_1(x -> 2*x) # Right-multiplying only. 
@test weights(2E_1) == 2*weights(E_1) # Left-multiplying only. 
@test 3E_1*z == (3E_1) * z == 3*(E_1 * z) # Linearity. 

# Test for error handling. 
testDist2 = Poisson(3)
@test_throws MethodError expectation(testDist2) # Catches unbounded distributions. 
@test_throws Exception E_1(x -> dot(x, ones(7))) # Test for non-applicable functions. MethodError on 0.6, DimensionMismatch on 0.7. 
@test_throws DimensionMismatch h.(z) * nodes(E_1) # Test non-commutativity of operator. 

#= 
    Tests for continuous distributions (no nodes supplied.)
=#

# Test Gauss-Legendre catchall for truncated normal. 
testDist3 = Truncated(Normal(0.5), 0, 10)
E_3 = expectation(testDist3)
mean3 = mean(testDist3)
@test E_3(x -> x) ≈ mean3 # First moment. 
@test E_3(x -> (x - mean3)^2) ≈ var(testDist3) # Second moment. 

# Test Gauss-Legendre catchall for arcsine distribution. (crazy)
glTestDist = Arcsine()
E_gl = expectation(glTestDist; n = 10^7) # Necessary because arcsine is crazy.
mean_gl = mean(glTestDist)
@test E_gl(x -> x) ≈ mean_gl atol = 1e-7
@test E_gl(x -> (x-mean_gl)^2) ≈ var(glTestDist) atol = 1e-7

# Test Gauss-Legendre catchall for triangular distribution. 
triTestDist = TriangularDist(1, 10)
E_tri = expectation(triTestDist, n = 1000)
@test E_tri(x -> x) ≈ mean(triTestDist) atol = 1e-5 # Necessary because of the discontinuity in the derivative. 

# Test for normal. 
testDist4 = Normal()
E_4 = expectation(testDist4)
@test E_4(x -> x) ≈ 0.0 atol = 1e-12 # First moment. 
@test E_4(x -> x^2) ≈ 1.0 atol = 1e-12 # Second moment. 
@test E_4(x -> x^3) ≈ 0.0 atol = 1e-12 # Third moment. 

# Test for lognormal. 
testDist5 = LogNormal()
E_5 = expectation(testDist5)
mean5 = mean(testDist5)
var5 = var(testDist5)
skew5 = skewness(testDist5)
@test E_5(identity) ≈ mean5 atol = 1e-12  # First moment. 
@test E_5(x -> (x-mean5)^2) ≈ var5 atol = 1e-12  # Second moment. 
@test E_5(x -> ((x-mean5)/sqrt(var5))^3) ≈ skew5 atol = 1e-12  # Third moment. 

# Test for beta.
testDist6 = Beta()
E_6 = expectation(testDist6)
mean6 = mean(testDist6)
var6 = var(testDist6)
skew6 = skewness(testDist6)
@test E_6(identity) ≈ mean6 atol = 1e-12  # First moment. 
@test E_6(x -> (x-mean6)^2) ≈ var6 atol = 1e-12  # Second moment. 
@test E_6(x -> ((x-mean6)/sqrt(var6))^3) ≈ skew6 atol = 1e-12  # Third moment. 

# Test error handling 
errorDist1 = Beta(Inf) # Beta Gaussian
@test_throws MethodError expectation(errorDist1)
errorDist2 = LogNormal(Inf, Inf) # LogNormal Gaussian
@test_throws MethodError expectation(errorDist2)
errorDist3 = Normal(0, Inf) # Normal Gaussian
@test_throws MethodError expectation(errorDist3)
errorDist4 = Uniform(-Inf, Inf) # Degenerate uniform, Gauss-Legendre. 
@test_throws MethodError expectation(errorDist4)

# Test some more distributions
gammaDist = Gamma()
expDist = Exponential()

mean_exp = mean(expDist)
var_exp = var(expDist)
skew_exp = skewness(expDist)
mean_gamma = mean(gammaDist)
var_gamma = var(gammaDist)
skew_gamma = skewness(gammaDist)

E_exp = expectation(expDist)
E_gamma = expectation(Gamma())

@test E_exp(x -> x) ≈ mean_exp atol = 1e-12 
@test E_gamma(x -> x) ≈ mean_gamma atol = 1e-12 
@test E_exp(x -> (x-1.0)^2) ≈ var_exp atol = 1e-12 
@test E_gamma(x -> (x-1.0)^2) ≈ var_gamma atol = 1e-12 
@test E_exp(x -> ((x-mean_exp)/sqrt(var_exp))^3) ≈ skew_exp atol = 1e-12 
@test E_gamma(x -> ((x-mean_gamma)/sqrt(var_gamma))^3) ≈ skew_gamma atol = 1e-12 

#= 
    Tests for some continuous distributions (nodes supplied.) 
=#

# Test trapezoidal for truncated normal. (REGULAR)
testDist7 = Truncated(Normal(), -5, 5)
z = -5:0.1:5
E_7 = expectation(testDist7, z)
E_8 = expectation(testDist7, z, Trapezoidal)
@test nodes(E_7) == nodes(E_8) && weights(E_7) == weights(E_8)
@test E_7(x -> x) ≈ mean(testDist7) atol = 1e-13
@test E_7(x -> x^2) ≈ var(testDist7) atol = 1e-6

# Test irregular trapezoidal.
z = unique([linspace(-5, 1, 203)' linspace(1, 5, 200)'])
E_9 = expectation(testDist7, z)
E_10 = expectation(testDist7, z, Trapezoidal)
@test nodes(E_9) == nodes(E_10) && weights(E_9) == weights(E_10)
@test E_9(x -> x) ≈ mean(testDist7) atol = 1e-8
@test E_9(x -> x^2) ≈ var(testDist7) atol = 1e-5

# Test convenience functions. 
@test expectation(x -> x, testDist7, z) == E_9(x -> x) # Grid one. 
@test expectation(x -> x, testDist4) == E_4(x -> x) # Standard normal/Gaussian continuous. 
@test expectation(x -> x, testDist) == E_1(x -> x) # Discrete one. 

# Test qnw dist results. 
@test expectation(x -> x, testDist4, QuantileLinSpace) ≈ 0.0 atol = 1e-12 
@test expectation(x -> x, testDist6, QuantileLinSpace) ≈ mean6 atol = 1e-12

# @test expectation(x -> ((x-mean_gamma)/sqrt(var_gamma))^3, gammaDist, QuantileLinSpace) ≈ skew_gamma atol = 1e-10