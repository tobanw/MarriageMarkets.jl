# MarriageMarkets

[![Build Status](https://travis-ci.org/tobanw/MarriageMarkets.jl.svg?branch=master)](https://travis-ci.org/tobanw/MarriageMarkets.jl)

The `MarriageMarkets` package currently provides two marriage market models as Julia types:

- `StaticMatch`: computes the equilibrium of the static frictionless marriage market model from "Who Marries Whom and Why" (Choo & Siow, 2006).
- `SearchMatch`: computes the equilibrium of variants on the search and matching model from "Assortative Matching and Search" (Shimer & Smith, 2000) and the empirical extension in "Marriage Gains" (Goussé, 2014).

`SearchMatch` also allows for inflows of new singles as well as deaths.

## Installation

As described in the manual, to [install unregistered packages][unregistered], use `Pkg.clone()` with the repository url:

```julia
Pkg.clone("git@github.com:tobanw/MarriageMarkets.jl.git")
```

Julia version 0.5 or higher is required (install instructions [here][version]).

## Usage

As `SearchMatch` supports a number of model variants, there are specific constructors for the two main types:

* `SearchClosed`: closed-system where agents cycle between singlehood and marriage
* `SearchInflow`: steady-state population is determined by exogenous inflows and type-specific death rates

## Example

The simple example below solves a search model with inflows and death.
I use [Gadfly][gadfly] to plot the match probability conditional on meeting.

```julia
using MarriageMarkets
using Gadfly
using Distributions

λ, δ = 500.0, 0.05 # arrival rates of meetings and divorce shocks
r = 0.05 # discount rate
σ = 1 # variance of Normally distributed match-specific productivity shocks
n = 50 # number of types
Θ = Vector(linspace(0.1, 1.0, n)) # types
f(x,y) = x*y # marital production function

γ = ones(n) ./ n # uniform inflows
ψ = ones(n) # uniform death rates


mgmkt = SearchInflow(λ, δ, r, σ, Θ, Θ, γ, γ, ψ, ψ, f)

plot(z=mgmkt.α, Geom.contour, Guide.title("Match probability conditional on meeting"))
```

![Match function](https://cloud.githubusercontent.com/assets/667531/23773990/8d362e86-04ef-11e7-8eb0-8d29b2b7ae96.png)

The saddle shape indicates positive assortative matching, as expected, due to the supermodular production function `f(x,y) = x*y`.

## Testing

In a Julia session, run `Pkg.test("MarriageMarkets")`.


[unregistered]:http://docs.julialang.org/en/release-0.5/manual/packages/#installing-unregistered-packages
[version]:http://julialang.org/downloads/platform.html
[gadfly]:http://gadflyjl.org/stable/
