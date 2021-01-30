# Mixture models
distset = [
    MixtureModel([Uniform(), Normal(), LogNormal()]),
    MixtureModel([Uniform(1, 10), Gamma()]),
    MixtureModel([Uniform(1, 100), Normal(1, 1000)]),
    MixtureModel([DiscreteUniform(1, 10), DiscreteUniform(1, 10)])
]

for dist in distset
    println(dist)
    μ = mean(dist)
    σ = std(dist)
    # No convenience call.
    E = expectation(dist)
    @test E(x -> x) ≈ μ
    @test E(x -> x^2) - μ^2 ≈ σ^2
    # Convenience call.
    @test expectation(x -> x, dist) ≈ μ
    @test expectation(x -> x^2, dist) - μ^2 ≈ σ^2
    # Stress tests.
    # Many nodes.
    E2 = expectation(dist, n = 100)
    @test E2(x -> x) ≈ μ
    @test E2(x -> x^2) - μ^2 ≈ σ^2
end

# Addition of expectations
e1 = expectation(Normal())
e2 = expectation(Uniform())
e3 = expectation(Gamma())

@test ((e1 + e2) + e3)(identity) ≈ ((e1 + e2) + e3)(identity) ≈ (e1 + e2 + e3)(identity)
@test (e1 + e2 + e3)(x -> x^2) ≈ 3.333333333333342

# Linear operator behavior.

# Error handling.
