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
@test_throws MethodError E_1(x -> dot(x, ones(7))) # Test for non-applicable functions. 
@test_throws DimensionMismatch h.(z) * nodes(E_1) # Test non-commutativity of operator. 

#= 
    Tests for some continuous distributions (no nodes supplied.)
=#

# Test Gauss-Legendre catchall for truncated normal. 
testDist3 = Truncated(Normal(0.5), 0, 10)
E_3 = expectation(testDist3)
@test E_3(x -> x) ≈ mean(testDist3) # Mean from mean(d::TruncatedNormal) in Distributions. 
squares = [x^2 for x in rand(testDist3, 10000000)] # 10^7
@test E_3(x -> x^2) ≈ mean(squares) atol = 1e-3

# Test Gauss-Legendre catchall for arcsine distribution. 
glTestDist = Arcsine()
E_gl = expectation(glTestDist)
@test E_gl(x -> x) ≈ mean(glTestDist) atol = 1e-3
squares = [x^2 for x in rand(glTestDist, 10000000)]
@test E_gl(x -> x^2) ≈ mean(squares) atol = 1e-3

# Test for normal. 
testDist4 = Normal()
E_4 = expectation(testDist4)
@test E_4(x -> sin(x)^2) ≈ 0.4320336833225541 atol = 1e-3 # From LLN average on 1 million points.

# Test for lognormal. 
testDist5 = LogNormal()
E_5 = expectation(testDist5)
@test E_5(x -> x) ≈ mean(testDist5) # Mean from mean (d::LogNormal) in Distributions. 
samples = rand(testDist5, 10000000) # 10^7 samples. 
sines = map(x -> sin(x)^3, samples)
@test_broken E_5(x -> sin(x)^3) ≈ sines # Broken for this kind of function. LLN n = 10^7. 
squares = map(x -> x^2, samples)
@test_broken E_5(x -> x^2) ≈ mean(squares) atol = 1e-4 

# Test for beta.
testDist6 = Beta()
E_6 = expectation(testDist6)
@test E_6(x -> x) ≈ mean(testDist6) # Mean from mean(d::LogNormal) in Distributions. 
samples = rand(testDist6, 10000000) # 10^7 samples.
squares = map(x -> x^2, samples)
@test E_6(x -> x^2) ≈ mean(squares) atol = 1e-3 # Allow for LLN error. 

# Test error handling 
errorDist1 = Beta(Inf) # Beta Gaussian
@test_throws MethodError expectation(errorDist1)
errorDist2 = LogNormal(Inf, Inf) # LogNormal Gaussian
@test_throws MethodError expectation(errorDist2)
errorDist3 = Normal(0, Inf) # Normal Gaussian
@test_throws MethodError expectation(errorDist3)
errorDist4 = Uniform(-Inf, Inf) # Degenerate uniform, Gauss-Legendre. 
@test_throws MethodError expectation(errorDist4)

#= 
    Tests for some continuous distributions (nodes supplied.) 
=#

# Test trapezoidal for truncated normal. (REGULAR)
testDist7 = Truncated(Normal(), -5, 5)
z = -5:0.1:5
E_7 = expectation(testDist7, z)
E_8 = expectation(testDist7, z, Trapezoidal)
@test nodes(E_7) == nodes(E_8) && weights(E_7) == weights(E_8)
@test E_7(x -> x) ≈ mean(testDist7) atol = 1e-10
squares = [x^2 for x in rand(testDist7, 10^7)]
@test E_7(x -> x^2) ≈ mean(squares) atol = 1e-3

# Test irregular trapezoidal.
z = unique([linspace(-5, 0, 100)' linspace(0, 5, 200)'])
E_9 = expectation(testDist7, z)
E_10 = expectation(testDist7, z, Trapezoidal)
@test nodes(E_9) == nodes(E_10) && weights(E_9) == weights(E_10)
@test E_9(x -> x) ≈ mean(testDist7) atol = 1e-4
@test E_9(x -> x^2) ≈ mean(squares) atol = 1e-3
