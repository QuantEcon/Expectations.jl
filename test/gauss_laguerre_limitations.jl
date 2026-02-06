# This file demonstrates the limitations of Gauss-Laguerre quadrature for
# Exponential and Gamma distributions, and shows how to use QuantileRange
# as an alternative for functions with variation near the mean.

using Test
using Expectations
using Distributions

@testset "Gauss-Laguerre limitations for Exponential" begin
    θ = 2.0
    dist = Exponential(θ)
    
    # Test 1: Verify that basic moments work well with Gauss-Laguerre
    E_gauss = expectation(dist, Gaussian, n=32)
    @test E_gauss(x -> x) ≈ mean(dist) rtol=1e-10  # Mean works well
    @test E_gauss(x -> x^2) ≈ 2*θ^2 rtol=1e-10     # Second moment works well
    
    # Test 2: Show that step functions have poor accuracy with Gauss-Laguerre
    # E[1_{X < θ}] = P(X < mean) = 1 - exp(-1) ≈ 0.632
    step_func(x) = x < θ ? 1.0 : 0.0
    analytical = 1 - exp(-1)
    
    result_gauss = E_gauss(step_func)
    # Gauss-Laguerre has significant error (>10% relative error)
    @test abs(result_gauss - analytical) > 0.05  # Error is substantial
    
    # Test 3: Show that QuantileRange provides better accuracy for step functions
    E_quant = expectation(dist, QuantileRange, n=50)
    result_quant = E_quant(step_func)
    # QuantileRange has much smaller error
    @test abs(result_quant - analytical) < 0.02  # Much better accuracy
    @test abs(result_quant - analytical) < abs(result_gauss - analytical)  # Better than Gauss-Laguerre
    
    # Test 4: Show node distribution issue
    # Gauss-Laguerre places most nodes in the tail
    nodes_gauss = nodes(E_gauss)
    @test sum(nodes_gauss .< θ) < 5  # Less than ~15% of nodes below mean
    @test sum(nodes_gauss .> θ) > 27  # More than ~85% of nodes above mean
    
    # QuantileRange distributes nodes more evenly
    nodes_quant = nodes(E_quant)
    @test sum(nodes_quant .< θ) > 20  # More balanced distribution
end

@testset "Gauss-Laguerre limitations for Gamma" begin
    α, θ = 2.0, 3.0
    dist = Gamma(α, θ)
    μ = mean(dist)  # α * θ = 6.0
    
    # Test 1: Verify that basic moments work well
    E_gauss = expectation(dist, Gaussian, n=32)
    @test E_gauss(x -> x) ≈ μ rtol=1e-10
    
    # Test 2: Show that functions with variation near mean have issues
    indicator_func(x) = abs(x - μ) < 1.0 ? 1.0 : 0.0
    
    result_gauss = E_gauss(indicator_func)
    
    # QuantileRange should provide better accuracy
    E_quant = expectation(dist, QuantileRange, n=50)
    result_quant = E_quant(indicator_func)
    
    # Both should be positive (some mass near the mean)
    @test result_gauss > 0
    @test result_quant > 0
    
    # QuantileRange typically gives different (often better) results for such functions
    # We don't assert which is better here as it depends on the specific function,
    # but we document that users should try both methods
end
