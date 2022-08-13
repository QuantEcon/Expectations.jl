# Common univariate distributions.
distset = [
    Normal(1.23, 3.45),
    Exponential(2.12),
    Gamma(4.3),
    Gamma(2.2, 3.1,),
    Erlang(2, 3),
    Beta(2.123),
    LogNormal(4.5),
    DiscreteUniform(1, 10),
    Chisq(3),
    Uniform(1, 200),
    Binomial(101, 0.45),
    Categorical([0.2, 0.35, 0.15, 0.3])
]

for dist in distset
    println(dist)
    μ = mean(dist)
    σ = std(dist)
    # No convenience call.
    E = expectation(dist)
    @test E(x -> x) ≈ μ
    @test E(x -> x^2) - μ^2 ≈ σ^2
    @test E(x -> ((x - μ)/σ)^3) + 1. ≈ skewness(dist) + 1. # To avoid comparisons to 0.0 exactly.
    # Convenience call.
    @test expectation(x -> x, dist) ≈ μ
    @test expectation(x -> x^2, dist) - μ^2 ≈ σ^2
    @test expectation(x -> ((x - μ)/σ)^3, dist) + 1. ≈ skewness(dist) + 1.
    # Stress tests.
    # Many nodes.
    E2 = expectation(dist, n = 100)
    @test E2(x -> x) ≈ μ
    @test E2(x -> x^2) - μ^2 ≈ σ^2
    @test E2(x -> ((x - μ)/σ)^3) + 1. ≈ skewness(dist) + 1. # To avoid comparisons to 0.0 exactly.
end

# Linear operator behavior.
distset = [
    DiscreteUniform(1., 10.),
    Normal(1.45),
    Beta(2.1)
]

for dist in distset
    E = expectation(dist)
    h(x) = 2*x
    z = nodes(E)
    @test E * h.(z) ≈ E(x -> 2*x) # Right-multiplying.
    @test weights(2E) ≈ 2*weights(E) # Left-multiplying.
    @test 3E*z ≈ (3E) * z ≈ 3*(E * z) # Linearity.
end

# Error handling.
distset = [ # Noncompact dists
    LogNormal(Inf, Inf),
    Beta(Inf),
    Normal(0, Inf),
    Uniform(-Inf, Inf),
    Poisson(3)
]

for dist in distset
    @test_throws ArgumentError expectation(dist)
end

distset = [ # Compact dist
    Uniform(1,2),
    Arcsine(1,2),
    Beta(1,1),
    Truncated(LogNormal(1,1), 0., 10.),
    ]

for dist in distset
    E = expectation(dist)
    h(x) = 2*x
    z = nodes(E)
    @test E * h.(z) ≈ E(x -> 2*x) # Right-multiplying.
    @test weights(2E) ≈ 2*weights(E) # Left-multiplying.
    @test 3E*z ≈ (3E) * z ≈ 3*(E * z) # Linearity.

    nodeList = nodes(E);
    E = expectation(dist, nodeList)
    @test E * h.(z) ≈ E(x -> 2*x) # Right-multiplying.
    @test weights(2E) ≈ 2*weights(E) # Left-multiplying.
    @test 3E*z ≈ (3E) * z ≈ 3*(E * z) # Linearity.
end


# Other errors.
E = expectation(DiscreteUniform(1, 10))
@test_throws Exception E(x -> dot(x, ones(7))) # Non-applicable functions.
@test_throws MethodError (x -> 2*x).(nodes(E)) * E # Non-commutativity.

# Trapezoidal methods.
distset = [Beta()]

for dist in distset
    # Setup.
    x = support(dist)
    μ = mean(dist)
    σ = std(dist)
    # Regular grid.
    grid = range(minimum(x), stop = maximum(x), length = 100)
    E = expectation(dist, grid)
    @test E(x -> x) ≈ μ
    @test abs(E(x -> x^2) - μ^2 - σ^2) < 1e-4
    # Irregular grid.
    grid2 = unique([grid' range(minimum(x), stop = maximum(x), length = 77)'])
    E2 = expectation(dist, grid2)
    @test E2(x -> x) isa Number # no accuracy guarantees for the irregular grid
    @test abs(E2(x -> x^2) - μ^2 - σ^2) isa Number
    # Convenience method
    @test expectation(identity, dist, grid2) ≈ E2(x -> x)
end

# Quantile

distset = [
    Uniform(1., 2.)
]

for dist in distset
    E = expectation(dist, QuantileRange)
    h(x) = 2*x
    z = nodes(E)
    @test E * h.(z) ≈ E(x -> 2*x) # Right-multiplying.
    @test weights(2E) ≈ 2*weights(E) # Left-multiplying.
    @test 3E*z ≈ (3E) * z ≈ 3*(E * z) # Linearity.
end

## Truncated distrubtions
@test_throws ArgumentError E = expectation(Pareto())
# Mean of pareto is (α * θ / (α - 1)). If we bound a Pareto at a high number we should get (close to) the analytical mean
α = 5.0
θ = 1.0
righttrunc = 10000
E = expectation(truncated(Pareto(α,θ),nothing,righttrunc),n=1000) # Right truncated Pareto at 1000
@test E(x->x) ≈ α*θ/(α-1)
