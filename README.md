[![Travis status](https://travis-ci.org/econtoolkit/Expectations.jl.svg?branch=master)](https://travis-ci.org/econtoolkit/Expectations.jl)

# Expectations

Installation:
```julia
Pkg.clone("http://github.com/econtoolkit/Expectations.jl.git")
```

This is a package designed to simplify the process of taking expectations of functions of random variables. The package is currently at **Julia v0.6**, and will break due to dependencies in v0.7.

### Random Variables 

The underlying distributions are objects from `Distributions.jl` (currently `<:UnivariateDistribution`).

### Quadrature Algorithms

There are two types of algorithm, `QuadratureAlgorithm` for algorithms which pick their own nodes, and 
`ExplicitQuadratureAlgorithm` for ones where the user picks. Currently, the only concrete subtypes of 
the former are `Gaussian` and `FiniteDiscrete`, and `Trapezoidal` for the latter (see "Mathematical Details"
for the mathematical details). 

### Expectation Operator

The package produces an expectation operators, `E`, as follows:

```julia
dist = Normal()
E = expectation(dist; kwargs...)
```

We can use these operators on functions as follows:

```julia
E(x -> x; kwargs...)
```

Expectations fall into a type hierarchy `IterableExpectation{NT, WT} <: IterableExpectation <: Expectation`
Currently, only iterables are supported, and they are parametrized by `weights` and `nodes`, and can therefore be used as 
follows:

```julia
h(x) = x^2 
dot(E.weights, h.(E.nodes))
```

Or as a linear operator, when defined with nodes:

```julia
dist = Normal()
z = -10:0.2:10
h = (x -> x^2).(z)
E = expectation(dist, nodes = z; kwargs...)
E*h # is equal to dot(h, E.weights)
3E*h # is equal to 3(E*h) and (3E)*h
```

### Mathematical Details

For finite discrete distributions, we simply compute the precise expectation. For 
continuous distributions given without nodes, we use some form of Gaussian quadrature
(either Gauss-Legendre, or a distribution-specific form for common distributions). These
distribution-specific algorithms are derived from `QuantEcon.jl`, who implemented them as 
defined in Miranda and Fackler's `CompEcon` toolbox. 

The default for distributions with nodes (and currently the only supported algorithm) is 
trapezoidal quadrature. 

It is important to ensure that the particular quadrature scheme used is compatible with the 
given function, distribution, and support. See the tests for accuracy under some normal and
edge cases. 

