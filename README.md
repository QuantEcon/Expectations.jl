[![Travis status](https://travis-ci.org/econtoolkit/Expectations.jl.svg?branch=master)](https://travis-ci.org/econtoolkit/Expectations.jl)

[![Travis 0.7]()

# Expectations

Installation:
```julia
Pkg.add("Expectations")
```

This is a package designed to simplify the process of taking expectations of functions of random variables. The package is compatible with Julia v0.6 and Julia v0.7. 

### Random Variables 

The underlying distributions are objects from `Distributions.jl` (currently `<:UnivariateDistribution`).

### Quadrature Algorithms

We have `QuadratureAlgorithm` for algorithms which pick their own nodes, and 
`ExplicitQuadratureAlgorithm` for ones where the user picks. Currently, the only concrete subtypes of 
the former are `Gaussian` and `FiniteDiscrete`, and `Trapezoidal` for the latter.

### Expectation Operator

The package produces an expectation operators, `E`, as follows:

```julia
dist = Normal()
E = expectation(dist)
E(x -> x)

E_morenodes = expectation(dist; n = 50) # Be careful, as too many nodes can introduce floating-point errors from miniscule exponents. 
```

Expectations fall into a type hierarchy `IterableExpectation{NT, WT} <: IterableExpectation <: Expectation`
Currently, only iterables are supported, and they are parametrized by `weights` and `nodes`.

```julia
h(x) = x^2 
dot(weights(E), h.(nodes(E)))
```

Or as a linear operator, either when you give your own nodes (or apply a function `h.(nodes(E))`):

```julia
dist = Normal()
z = -10:0.2:10
h = (x -> x^2).(z)
E = expectation(dist, z; kwargs...)
E*h # is equal to dot(h, weights(E))
3E*h # is equal to 3(E*h) and (3E)*h
```

### Mathematical Details

For finite discrete distributions, we simply compute the precise expectation. For 
continuous distributions given without nodes, we use some form of Gaussian quadrature
(either Gauss-Legendre, or a distribution-specific form for common distributions). These
distribution-specific algorithms are derived from `FastGaussQuadrature.jl`.

The default for distributions with nodes (and currently the only supported algorithm) is 
trapezoidal quadrature. 

It is important to ensure that the particular quadrature scheme used is compatible with the 
given function, distribution, and support. See the tests for accuracy under some normal and
edge cases. 

