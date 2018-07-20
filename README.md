[![Travis status](https://travis-ci.org/econtoolkit/Expectations.jl.svg?branch=master)](https://travis-ci.org/econtoolkit/Expectations.jl)
[![codecov](https://codecov.io/gh/econtoolkit/Expectations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/econtoolkit/Expectations.jl)
[![Coverage Status](https://coveralls.io/repos/github/econtoolkit/Expectations.jl/badge.svg?branch=master)](https://coveralls.io/github/econtoolkit/Expectations.jl?branch=master)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://econtoolkit.github.io/Expectations.jl/latest)

# Expectations

Installation:
```julia
Pkg.add("Expectations")
```

This is a package designed to simplify the process of taking expectations of functions of random variables. The package is compatible with Julia v0.6 and Julia v0.7. 

### Random Variables 

The underlying distributions are objects from `Distributions.jl` (currently `<:UnivariateDistribution`).

### Quadrature Algorithms

We support different types of Gaussian quadrature (Gauss-Hermite, Gauss-Legendre, Gauss-Laguerre, etc.) based on the distribution, as well as some methods
with user-defined nodes (e.g., trapezoidal integration).

### Expectation Operator

The key object is the expectation operator, `E`, which can be used as follows:

```julia
dist = Normal()
E = expectation(dist)
E(x -> x)
```

Or as a linear operator on vectors: 

```julia
dist = Normal()
z = -10:0.2:10
h = (x -> x^2).(z)
E = expectation(dist, z; kwargs...)
E*h # is equal to dot(h, weights(E))
3E*h # is equal to 3(E*h) and (3E)*h
```
