using Distributions

ρ = 5.0
r = 0.05
δ = 0.05
σ = 10.0

const STDNORMAL = Normal()
h(x::Real, y::Real) = x*y # unidimensional production function
hm(x,y) = dot(x,y)/2 # multidimensional production function
