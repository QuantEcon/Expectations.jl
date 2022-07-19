## Overview

The goal of this package is to provide an intuitive and mathematically sound interface for taking expectations of random variables
and their higher-order functions (i.e., if ``X \sim N(0, 1)``, what is ``\mathbb{E}[\sin(X)]``?).

The underlying random variables are assumed to be distributions from [`Distributions.jl`](https://github.com/juliastats/distributions.jl). Currently, only univariate distributions are supported.

## Installation

To install, run:

```@repl 1
using Pkg
pkg"add Expectations Distributions"

using Expectations, Distributions
```

## The Expectation Operator

The key object in this package is an **expectation operator**, or an object `<: Expectation`. These include all objects capable of being called on a function; e.g. that support a method `function (e::Expectation)(f::Function)`. You can create these as following:

```@repl 1

dist = Normal();
E = expectation(dist)
```

You can also choose algorithms and default parameters (see below for list):

```@repl 1
E = expectation(dist, Gaussian; n = 30) # Could have done expectation(dist) or expectation(dist; n = 30)
```

These objects can then be applied to functions:

```@repl 1
E(x -> x)
E(x -> x^2)
```

There is also a convenience function to evaluate expectations directly, without returning the operator:

```@repl 1
f = x -> x^2
expectation(f, dist)
```

In general, `expectation(f, dist, ...)` is equivalent to `E(f)`, where `E = expectation(dist, ...)`.

### IterableExpectation

The only concrete subtype of `Expectation` currently supported is `IterableExpectation{NT, WT}`. These are expectations for which we have a
discrete vector of quadrature nodes and weights, either defined by user fiat, or set algorithmically. These support some additional behavior:

```@repl 1
nodeList = nodes(E);
vals = map(x -> x^2, nodeList);
E * vals
(2E) * vals
```

The above behavior, in some sense, puts the "operator" in "expectation operator"; that is, it allows it to move elements of a vector space around, and to be scalar-multiplied.

### User-Defined Nodes

There are some situations where we are forced to use a specific set of nodes. In those situations, `E = expectation(dist, nodes)` will create the relevant object.

## Mixture Models

We also have support for univariate [mixture models](https://juliastats.org/Distributions.jl/latest/mixture). The `MixtureExpectation` type is a struct with two fields:

* `expectations`, which is a list of `IterableExpectation` objects.

* `mixtureweights`, which is the mixing probabilities over the various model components.

The mixture models are constructed using the **default settings** for each component distribution. (It still accepts kwargs, which are applied to each.)

```@repl 1
d = MixtureModel([Uniform(), Normal(), Gamma()]);
E = expectation(d; n = 30); # n = 30 nodes for each
@show typeof(E)
E(x -> abs(x))
```

If you want to change this, you should construct each distribution separately, and then chain them together.

```@repl 1
E1 = expectation(Uniform())
E2 = expectation(Normal())
E3 = expectation(Gamma())

E = MixtureExpectation([E1, E2, E3], [1/3, 1/3, 1/3])
E(x -> abs(x))
```

## Supported Distributions, Algorithms, Keywords, and Defaults

Here is a list of currently supported distributions, along with keyword arguments and their defaults.  

| Distribution Name | Algorithm (Julia Type) | Keywords and Defaults | Restrictions |
| ----------------- | -------------- | --------------------- | ------------ |
| Discrete Univariate | FiniteDiscrete <: QuadratureAlgorithm | N/A | Support must be finite. |
| Continuous Univariate | Gauss-Legendre (Gaussian <: QuadratureAlgorithm) | n = 500 | Support must be a compact interval ``[a, b]``. |
| Continuous Univariate | QNWDist[^1] (QuantileRange <: ...) | n = 50, q0 = 0.001, qN = 0.999 | Distribution must be nondegenerate. |
| Normal <: Continuous Univariate | Gauss-Hermite (...) | n = 30 | ... |
| LogNormal <: ... | Gauss-Hermite (...) | n = 30 | ... |
| Beta <: ... | Gauss-Jacobi (...) | n = 32 | ... |
| ChiSq <: ... | Gauss-Laguerre (...) | n = 32 | ... |
| Uniform <: ... | Gauss-Legendre (...) | n = 30 | ... |
| Exponential <: ... | Gauss-Laguerre (...) | n = 32 | ... |
| Gamma <: ... | Gauss-Laguerre (...) | n = 32 | ... |
| Univariate | Trapezoidal <: ExplicitQuadratureAlgorithm | N/A | All nodes must be inside distribution's support. |

###

Some unbounded distributions are currently not supported (e.g., Poisson). Depending on your use case truncating may be a feasible option:
```@repl
E = expectation(Pareto()) # Throws error
E = expectation(truncated(Pareto(),0.0,1000.0)) # Truncated Pareto on [0,1000]
```
See `Distributions.truncated` for more. Of course, truncating the distribution also affects its properties.

## Mathematical Details and References

The specific quadrature algorithms come from the [`FastGaussQuadrature.jl`](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) library, which is maintained by [Alex Townsend](https://github.com/ajt60gaibb) of Cornell University. Much of the quadrature code came from the [`DistQuads.jl`](https://github.com/pkofod/DistQuads.jl) library, which is maintained by [Patrick K. Mogensen](https://github.com/pkofod) at the University of Copenhagen. In addition, there are some objects contributed by individual users; see docstring for citations.

> **WARNING**: It is important to be aware of the deficiencies of numerical quadrature schemes. For example, it is recommended to be careful when using these methods for the following classes of functions and situations:

* Discontinuous or nondifferentiable functions (even if the function is a.e.-differentiable)
* Periodic/oscillatory functions with a high frequency
* Extremely large numbers of quadrature nodes, which may lead to vanishingly small weights.

## Contact

If you would like to get in touch, please do one of the following:

* Issue requests: Open an issue on the [package repository](https://github.com/QuantEcon/Expectations.jl) with the tag `feature request`.
* Bugs: Same as above, but with the tag `bug`.
* Pull Request: We are always open to new functionality. If you have a feature you'd like to add (say, a new distribution or algorithm), once you prepare a PR with the feature and some tests, open it in the usual way.
* Other: You can reach out to Arnav Sood at [`misc@arnavsood.com`](mailto:misc@arnavsood.com)
* Citation: If this package was helpful in your research work, you may consider citing the package in whatever method is appropriate for your field.

[^1]: This is a quadrature scheme written by [Spencer Lyon](http://spencerlyon.com/) (PhD. NYU) as part of the [`QuantEcon`](https://quantecon.org/) project. Used with permission.
