[![Travis status](https://travis-ci.org/econtoolkit/Expectations.jl.svg?branch=master)](https://travis-ci.org/econtoolkit/Expectations.jl)
[![codecov](https://codecov.io/gh/econtoolkit/Expectations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/econtoolkit/Expectations.jl)
[![Coverage Status](https://coveralls.io/repos/github/econtoolkit/Expectations.jl/badge.svg?branch=master)](https://coveralls.io/github/econtoolkit/Expectations.jl?branch=master)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://econtoolkit.github.io/Expectations.jl/latest)

# Expectations

Installation (v0.7 and up):
```julia
pkg> add Expectations
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
For convenience,
```julia
expectation(x->x^2, dist)
```

As a linear operator on vectors using the nodes of the distribution 
```julia
dist = Normal()
E = expectation(dist)
x = nodes(E)
f(x) = x^2
E * f.(x) == dot(f.(x), weights(E))
```

If nodes are given, it will calculate using Newton-Coates quadrature (e.g. Trapezoidal)
```julia
x = -10:0.2:10
f(x) = x^2
E = expectation(dist, x)
3 * E(f) == 3 * E * f.(x)
```
