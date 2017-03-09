### Verifications: basic Shimer-Smith model ###

"""
Element-by-element calculation of steady state flow equations

```julia
mres = (δ .+ ψ_m) .* ℓ_m - u_m .* ((δ .+ ψ_m) + ρ * (α * u_f))
fres = (δ .+ ψ_f) .* ℓ_f - u_f .* ((δ .+ ψ_f) + ρ * (α' * u_m))
```
"""
function sse_base(M::SearchMatch)

	mres = similar(M.ℓ_m)
	fres = similar(M.ℓ_f)
	
	for i in 1:length(M.ℓ_m)
		mres[i] = (M.δ + M.ψ_m[i]) * M.ℓ_m[i] - M.u_m[i] * ((M.δ + M.ψ_m[i]) + M.ρ * dot(M.α[i,:], M.u_f))
	end
	for j in 1:length(M.ℓ_f)
		fres[j] = (M.δ + M.ψ_f[j]) * M.ℓ_f[j] - M.u_f[j] * ((M.δ + M.ψ_f[j]) + M.ρ * dot(M.α[:,j], M.u_m))
	end
	return mres, fres
end

"""
Element-by-element calculation of value function equations
```julia
mres = 2*v_m - ρ * (αS * u_f)
fres = 2*v_f - ρ * (αS' * u_m),
```
where `αS = α .* s/(r+δ+ψ_m(x)+ψ_f(y))`.
"""
function vf_base(M::SearchMatch)

	mres = similar(M.v_m)
	fres = similar(M.v_f)

	for i in 1:length(M.v_m)
		mres[i] = 2*M.v_m[i] - M.ρ * dot(M.α[i,:] .* M.s[i,:] ./ (M.r+M.δ+M.ψ_m[i].+M.ψ_f), M.u_f)
	end
	for j in 1:length(M.v_f)
		fres[j] = 2*M.v_f[j] - M.ρ * dot(M.α[:,j] .* M.s[:,j] ./ (M.r+M.δ+M.ψ_f[j].+M.ψ_m), M.u_m)
	end
	return mres, fres
end

function surplus_base(M::SearchMatch)
	s = similar(M.h)

	for i in 1:length(M.v_m), j in 1:length(M.v_f)
		s[i,j] = M.h[i,j] - M.v_m[i] - M.v_f[j]
	end

	return s
end


### Setup ###

ρ = 5.0
δ = 0.05
r = 0.05
σ = 0

# symmetric case only
ntypes = 50
mass = 100
types = Vector(linspace(1.0, 2.0, ntypes))

# uniform population distributions
lm = (mass / ntypes) * ones(Float64, ntypes)
lf = (mass / ntypes) * ones(Float64, ntypes)

# production function
hsup(x::Real, y::Real) = x*y

# instantiate a MarriageMarket
symm = SearchClosed(ρ, δ, r, σ, types, types, lm, lf, hsup)


### Tests: basic Shimer-Smith model ###

# check convergence
msse, fsse = sse_base(symm)
mvf, fvf = vf_base(symm)
α_err = symm.α .- convert(Array{Float64}, (surplus_base(symm) .> 0.0))

# valid solution
@fact msse --> roughly(zeros(msse), atol = 1e-7)
@fact fsse --> roughly(zeros(fsse), atol = 1e-7)
@fact mvf --> roughly(zeros(mvf), atol = 1e-7)
@fact fvf --> roughly(zeros(fvf), atol = 1e-7)
@fact α_err --> zeros(α_err)

# symmetry
@fact symm.v_m --> roughly(symm.v_f)
@fact symm.u_m --> roughly(symm.u_f)
@fact symm.α --> roughly(symm.α')
