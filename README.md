# MarriageMarkets

[![Build Status](https://travis-ci.org/tobanw/MarriageMarkets.jl.svg?branch=master)](https://travis-ci.org/tobanw/MarriageMarkets.jl)

The `MarriageMarkets` package currently provides two marriage market models as Julia types:

- `ChooSiow`: computes the equilibrium of the static frictionless marriage market model from "Who Marries Whom and Why," Choo & Siow (2006).
- `ShimerSmith`: computes the equilibrium of the search and matching model from "Assortative Matching and Search," Shimer & Smith (2000).

`ShimerSmith` is generalized in a few ways beyond the model presented in the paper.
This includes a match-specific "love" shock to make the matching probabilistic.

Upcoming features:

- exogenous inflow of singles and outflow via death
- type depreciation (to capture aging, for example)
- multi-dimensional types

## Installation

As described in the manual, to [install unregistered packages][install], use `Pkg.clone()` with the repository url:

```julia
Pkg.clone("git@github.com:tobanw/MarriageMarkets.jl.git")
```

## Usage

Look at the unit tests for examples of using the types provided by `MarriageMarkets`.

## Testing

In a Julia session, run `Pkg.test("MarriageMarkets")`.


[install]: http://docs.julialang.org/en/release-0.5/manual/packages/#installing-unregistered-packages
