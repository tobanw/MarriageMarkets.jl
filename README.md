# MarriageMarkets

[![Build Status](https://travis-ci.org/tobanw/MarriageMarkets.jl.svg?branch=master)](https://travis-ci.org/tobanw/MarriageMarkets.jl)

The `MarriageMarkets` package currently provides two marriage market models as Julia types:

- `StaticMatch`: computes the equilibrium of the static frictionless marriage market model from "Who Marries Whom and Why," Choo & Siow (2006).
- `SearchMatch`: computes the equilibrium of variants on the search and matching model from "Assortative Matching and Search," Shimer & Smith (2000).

`SearchMatch` is generalized in a few ways beyond the model presented in the paper.
This includes:

- a match-specific "love" shock to make the matching probabilistic
- exogenous inflow of singles and outflow via death

Upcoming features:

- type depreciation (to capture aging, for example)
- endogenous divorce
- multi-dimensional types

## Installation

As described in the manual, to [install unregistered packages][unregistered], use `Pkg.clone()` with the repository url:

```julia
Pkg.clone("git@github.com:tobanw/MarriageMarkets.jl.git")
```

Julia version 0.5 or higher is required (install instructions [here][version]).

## Usage

As `SearchMatch` supports a number of model variants, there are some convenience methods for the main model specifications.

* `SearchClosed`: closed-system where agents cycle between singlehood and marriage
* `SearchInflow`: steady-state population is determined by exogenous inflows and type-specific death rates

Look at the unit tests for more examples of using the types provided by `MarriageMarkets`.

## Example

The simple example below solves a search model with inflows and death.
I use [Gadfly][gadfly] to plot the match probability conditional on meeting.

```julia
using MarriageMarkets
using Gadfly
using Distributions

n = 50 # number of types
ρ, δ, r = 500, 0.05, 0.05
Θ = Vector(linspace(0.0, 1.0, n)) # types
h(x,y) = x*y # marital production function

γ = ones(n) ./ n # uniform inflows
ψ = ones(n) # uniform death rates

G(x) = cdf(Normal(0.0, 0.1), x) # distribution of shocks

mgmkt = SearchInflow(ρ, δ, r, Θ, Θ, γ, γ, ψ, ψ, h, G)

plot(z=mgmkt.α, Geom.contour)
```

![Match function](https://cloud.githubusercontent.com/assets/667531/20243457/fca3c682-a925-11e6-9c9e-2baee549d7bb.png)

## Testing

In a Julia session, run `Pkg.test("MarriageMarkets")`.


[unregistered]:http://docs.julialang.org/en/release-0.5/manual/packages/#installing-unregistered-packages
[version]:http://julialang.org/downloads/platform.html
[gadfly]:http://gadflyjl.org/stable/
