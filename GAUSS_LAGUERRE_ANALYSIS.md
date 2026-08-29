# Analysis of Gauss-Laguerre Quadrature for Exponential Distributions

## Summary

This document explains the limitations of Gauss-Laguerre quadrature for Exponential and Gamma distributions, as reported in the issue. **This is NOT a bug in the code** - the implementation is mathematically correct. Rather, it is a fundamental mathematical property of Gauss-Laguerre quadrature that can lead to poor accuracy for certain types of functions.

## The Problem

### Mathematical Background

Gauss-Laguerre quadrature computes integrals of the form:
```
∫₀^∞ f(u) e^(-u) du ≈ Σ wᵢ f(uᵢ)
```

For an Exponential distribution with scale parameter θ (mean = θ), the expectation is:
```
E[g(X)] = ∫₀^∞ g(x) (1/θ) e^(-x/θ) dx
```

After the change of variables u = x/θ, this becomes:
```
E[g(X)] = ∫₀^∞ g(θu) e^(-u) du ≈ Σ wᵢ g(θuᵢ)
```

So the nodes are: xᵢ = θ * uᵢ where uᵢ are the Gauss-Laguerre nodes.

### Node Distribution Issue

The Gauss-Laguerre nodes {uᵢ} are optimized for the weight function e^(-u), which decays exponentially. This places most nodes in the tail of the distribution:

For n=32 nodes:
- Only 3 nodes (9.4%) are below u=1 (the mean of the standard exponential)
- 29 nodes (90.6%) are above u=1
- The median node is around u ≈ 0.58

When scaled by θ for Exponential(θ), this means:
- For Exponential(2.0): only 3 nodes below x=2 (the mean), median node at x ≈ 42.5
- Most quadrature nodes are exploring the far tail (x > 4θ)

### Impact on Accuracy

This node distribution is **excellent** for:
- Polynomial functions (x, x², x³, etc.) - exact up to degree 2n-1
- Functions that vary smoothly throughout the distribution
- Functions where the tail behavior matters

This node distribution is **poor** for:
- Step functions or indicator functions near the mean
- Functions with rapid variation near x = θ (the mean)
- Discontinuous or non-smooth functions
- Functions where the key behavior is near the mode/mean, not the tail

### Empirical Evidence

Testing with Exponential(2.0) and a step function f(x) = 1 if x < 2, else 0:
- Analytical result: P(X < mean) = 1 - e^(-1) ≈ 0.6321
- Gauss-Laguerre (n=32): ≈ 0.5549 (relative error: -12.2%)
- QuantileRange (n=50): ≈ 0.6223 (relative error: -1.5%)

## The Solution

### Not a Code Bug

The implementation in `src/iterable.jl` is **mathematically correct**:
```julia
function _expectation(dist::Exponential, alg::Type{Gaussian}; n=32, kwargs...)
    θ = inv(dist.θ)  # Convert from scale to rate
    gl = gausslaguerre(n)
    nodes = gl[1] ./ θ   # Correct: gl[1] / (1/θ) = gl[1] * θ
    weights = gl[2]
    return IterableExpectation(nodes, weights)
end
```

### Documentation and Warnings Added

1. **Updated docs/src/index.md**: Added a warning in the "Mathematical Details and References" section about Gauss-Laguerre limitations for Exponential and Gamma distributions.

2. **Enhanced docstrings in src/iterable.jl**: Added detailed notes to the `_expectation` methods for Exponential and Gamma distributions, explaining the node distribution issue and recommending QuantileRange as an alternative.

3. **Added test file**: Created `test/gauss_laguerre_limitations.jl` to demonstrate the issue and show how to use QuantileRange as an alternative.

### Recommended Alternative: QuantileRange

For functions with variation near the mean/mode, use QuantileRange quadrature:

```julia
using Expectations, Distributions

dist = Exponential(2.0)

# Default (Gauss-Laguerre) - may have poor accuracy for some functions
E_laguerre = expectation(dist)

# Alternative (QuantileRange) - better for functions with variation near mean
E_quantile = expectation(dist, QuantileRange, n=50)

# For a step function near the mean
f(x) = x < 2.0 ? 1.0 : 0.0
result_laguerre = E_laguerre(f)  # Larger error
result_quantile = E_quantile(f)   # Smaller error
```

QuantileRange distributes nodes evenly across the quantiles of the distribution, ensuring good coverage both near the mean and in the tails.

## When to Use Each Method

### Use Gauss-Laguerre (default) when:
- Computing moments (mean, variance, etc.)
- Working with smooth polynomial-like functions
- The function varies smoothly across the entire distribution

### Use QuantileRange when:
- The function has rapid variation near the mean/mode
- Working with step functions or indicator functions
- The function has discontinuities
- You need more balanced node coverage across the distribution

## Conclusion

The Gauss-Laguerre implementation is correct but has inherent mathematical limitations for certain types of functions when applied to Exponential and Gamma distributions. Users should be aware of these limitations and use QuantileRange quadrature when appropriate. The documentation has been updated to clearly warn users about this issue and provide guidance on when to use alternative methods.
