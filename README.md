# MarriageMarkets

[![Build Status](https://travis-ci.org/tobanw/MarriageMarkets.jl.svg?branch=master)](https://travis-ci.org/tobanw/MarriageMarkets.jl)

The `MarriageMarkets` package currently provides two marriage market models as Julia types.
First is the `ChooSiow` type, which computes the equilibrium of the static frictionless marriage market model from "Who Marries Whom and Why," Choo & Siow (2006).
Second is the `ShimerSmith` type, which computes the equilibrium of the dynamic search and matching model from "Assortative Matching and Search," Shimer & Smith (2000).

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
