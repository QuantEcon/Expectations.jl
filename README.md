[![Travis status](https://travis-ci.org/QuantEcon/Expectations.jl.svg?branch=master)](https://travis-ci.org/QuantEcon/Expectations.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/i67ucrxj4yf1kdx8?svg=true)](https://ci.appveyor.com/project/arnavs/expectations-jl)
[![codecov](https://codecov.io/gh/QuantEcon/Expectations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantEcon/Expectations.jl)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://QuantEcon.github.io/Expectations.jl/dev)

# Expectations

Installation (v1.0 and up):
```julia
pkg> add Expectations
```

This is a package designed to simplify the process of taking expectations of functions of random variables. 

### Expectation Operator

The key object is the `expectation` function, which returns an operator:

```julia
dist = Normal()
E = expectation(dist)
E(x -> x)
```
For convenience, the operator can be applied directly to a function instead of being cached,
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

### Random Variables 

The underlying distributions are objects from `Distributions.jl` (currently `<:UnivariateDistribution`).

### Quadrature Algorithms

We support different types of Gaussian quadrature (Gauss-Hermite, Gauss-Legendre, Gauss-Laguerre, etc.) based on the distribution, as well as some methods with user-defined nodes (e.g., trapezoidal integration).

We have rules for the following distributions: 

* Normal 
* ChiSq
* LogNormal
* Exponential
* Beta
* Gamma/Erlang
* Continuous Univariate (compact; generic fallback)
* Continuous Univariate (no restriction; approximates with quantile grid)
* Discrete

See docs for more info. 
