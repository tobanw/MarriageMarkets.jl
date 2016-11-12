"""
A Julia package for solving various marriage market matching models.

Two models are supported:

- `ChooSiow`: the static frictionless matching model of Choo & Siow (2006)
- `ShimerSmith`: the search and matching model of Shimer & Smith (2000)

`ShimerSmith` is generalized in a few ways beyond the model presented in the paper.
This includes a match-specific "love" shock to make the matching probabilistic.

Upcoming features:

- exogenous inflow of singles and outflow via death
- type depreciation (to capture aging, for example)
- multi-dimensional types
"""
module MarriageMarkets

export ChooSiow, ShimerSmith, SearchClosed, SearchInflow

include("choo-siow.jl")
include("shimer-smith.jl")

end # module
