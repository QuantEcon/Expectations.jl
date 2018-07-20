## Overview 

The goal of this package is to provide an intuitive and mathematically sound interface for taking expectations of random variables
and their higher-order functions (i.e., if ``X \sim N(0, 1)``, what is ``\mathbb{E}[\sin(X)]``?). 

The underlying random variables are assumed to be distributions from [`Distributions.jl`](https://github.com/juliastats/distributions.jl). Currently, 
only univariate distributions are supported. 

## Installation 

To install, run (in 0.7): 

```@repl 
using Expectations
```

Currently, Julia 0.6 and Julia 0.7 are supported. 

## The Expectation Operator 

The key object in this package is an **expectation operator**, or an object `<: Expectation`. These include all objects capable of being called on a function; e.g. that support a method `function (e::Expectation)(f::Function)`. You can create these as following:

```@repl
dist = Normal();
E = expectation(dist)
```

You can also choose and algorithms and default parameters (see below for list):

```@repl
E = expectation(dist, Gaussian; n = 30)
```

These objects can then be applied to functions: 

```@repl
E(x -> x)
E(x -> x^2)
```

### IterableExpectation

The only concrete subtype of `Expectation` currently supported is `IterableExpectation{NT, WT}`. These are expectations for which we have a
discrete vector of quadrature nodes and weights, either defined by user fiat, or set algorithmically. These support some additional behavior: 

```@repl
nodeList = nodes(E)
vals = map(x -> sin(x)^2, nodeList)
E * vals
(2E) * vals
```

The above behavior, in some sense, puts the "operator" in "expectation operator"; that is, it allows it to move elements of a vector space around, 
and to be scalar-multiplied. 


## Supported Distributions, Algorithms, Keywords, and Defaults 

Here is a list of currently supported distributions, along with keyword arguments and their defaults.  

| Distribution Name | Algorithm (Julia Type) | Keywords and Defaults | Restrictions | 
| ----------------- | -------------- | --------------------- | ------------ | 
| Discrete Univariate | FiniteDiscrete <: QuadratureAlgorithm | N/A | Support must be finite. | 
| Continuous Univariate | Gauss-Legendre (Gaussian <: QuadratureAlgorithm) | n = 500 | Support must be a compact interval ``[a, b]``. |
| Normal <: Continuous Univariate | Gauss-Hermite (...) | n = 30 | Parameters must be finite. | 
| LogNormal <: ... | Gauss-Hermite (...) | n = 30 | ... | 
| Beta <: ... | Gauss-Jacobi (...) | n = 32 | ... | 
| Exponential <: ... | Gauss-Laguerre (...) | n = 32 | ... | 
| Gamma <: ... | Gauss-Laguerre (...) | n = 32 | ... | 
| Univariate | Trapezoidal <: ExplicitQuadratureAlgorithm | N/A | All nodes must be inside distribution's support. | 

## Mathematical Details and References 

The specific quadrature algorithms come from the [`FastGaussQuadrature.jl`](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) library, which is maintained by [Alex Townsend](https://github.com/ajt60gaibb) of Cornell University. Much of the quadrature code came from the [`DistQuads.jl`](https://github.com/pkofod/DistQuads.jl) library, which is maintained by [Patrick K. Mogensen](https://github.com/pkofod) at the University of Copenhagen. 

:warning: It is important to be aware of the deficiencies of numerical quadrature schemes. For example, it is recommended to be careful when using these methods for the following classes of functions and situations: 

* Discontinuous or nondifferentiable functions (even if the function is a.e.-differentiable)
* Periodic/oscillatory functions with a high frequency 
* Extremely large numbers of quadrature nodes, which may lead to vanishingly small weights. 

## Contact 

If you would like to get in touch, please do one of the following:

* Issue requests: Open an issue on the [package repository](https://github.com/econtoolkit/Expectations.jl) with the tag `feature request`. 
* Bugs: Same as above, but with the tag `bug`. 
* Pull Request: We are always open to new functionality. If you have a feature you'd like to add (say, a new distribution or algorithm), once you prepare a PR with the feature and some tests, open it in the usual way. 
* Other: You can reach out to Jesse Perla at [`jesse.perla@ubc.ca`](mailto:jesse.perla@ubc.ca) and Arnav Sood at [`arnav.sood@ubc.ca`](mailto:arnav.sood@ubc.ca)
* Citation: If this package was helpful in your research work, you may consider citing the package in whatever method is appropriate for your field. 