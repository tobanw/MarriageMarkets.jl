ρ = 5.0
r = 0.05
δ = 0.05

function closed_uniform(ntypes, mmass, fmass, prod, σ)
	
	# types
	Θ = Vector(linspace(1.0, 2.0, ntypes))

	# uniform population distributions
	lm = (mmass / ntypes) * ones(Float64, ntypes)
	lf = (fmass / ntypes) * ones(Float64, ntypes)

	return SearchClosed(ρ, δ, r, σ, Θ, Θ, lm, lf, prod)
end

function inflow_uniform(ntypes, mmass, fmass, prod, σ)
	
	# types
	Θ = Vector(linspace(1.0, 2.0, ntypes))

	# uniform population inflows
	gm = (mmass / ntypes) * ones(Float64, ntypes)
	gf = (fmass / ntypes) * ones(Float64, ntypes)

	# uniform death rates: normalized to 1
	d = ones(Float64, ntypes)

	return SearchInflow(ρ, δ, r, σ, Θ, Θ, gm, gf, d, d, prod)
end

h(x::Real, y::Real) = x*y

"""
Element-by-element calculation of

```julia
mres = (M.δ .+ M.ψ_m) .* M.ℓ_m - M.u_m .* ((M.δ .+ M.ψ_m) + M.ρ * (M.α * M.u_f))
fres = (M.δ .+ M.ψ_f) * M.ℓ_f - M.u_f .* ((M.δ .+ M.ψ_f) + M.ρ * (M.α' * M.u_m))
```
"""
function sse_resid(M::SearchMatch)

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
Element-by-element calculation of
		``∀x, 2w(x) = ρ ∫ α(x,y) S(x,y) u(y) dy``, where

```julia
mres = 2*M.w_m - ρ * (αS * M.u_f)
fres = 2*M.w_f - ρ * (αS' * M.u_m),
```
where `αS = M.α .* M.S`.
"""
function vf_resid(M::SearchMatch)

	mres = similar(M.w_m)
	fres = similar(M.w_f)

	for i in 1:length(M.w_m)
		mres[i] = 2*M.w_m[i] - M.ρ * dot(M.α[i,:] .* M.S[i,:], M.u_f)
	end
	for j in 1:length(M.w_f)
		fres[j] = 2*M.w_f[j] - M.ρ * dot(M.α[:,j] .* M.S[:,j], M.u_m)
	end
	return mres, fres
end

function match_surplus(M::SearchMatch)
	S = similar(M.h)

	for i in 1:length(M.w_m), j in 1:length(M.w_f)
		S[i,j] = (M.h[i,j] - M.w_m[i] - M.w_f[j]) / (M.r + M.δ + M.ψ_m[i] + M.ψ_f[j])
	end

	return S
end


### Deterministic case

# symmetric case only
symm = closed_uniform(50, 100, 100, h, 0)

# check convergence
msse, fsse = sse_resid(symm)
mvf, fvf = vf_resid(symm)
α_err = symm.α .- convert(Array{Float64}, (match_surplus(symm) .> 0.0))

# valid solution
@fact msse --> roughly(zeros(msse), atol = 1e-7)
@fact fsse --> roughly(zeros(fsse), atol = 1e-7)
@fact mvf --> roughly(zeros(mvf), atol = 1e-7)
@fact fvf --> roughly(zeros(fvf), atol = 1e-7)
@fact α_err --> zeros(α_err)

# symmetry
@fact symm.w_m --> roughly(symm.w_f)
@fact symm.u_m --> roughly(symm.u_f)
@fact symm.α --> roughly(symm.α')


### Randomness: marital productivity shocks

using Distributions

σ = 10
G(x) = cdf(Normal(0, σ), x)

# symmetric case
rsym = closed_uniform(50, 100, 100, h, σ)

@fact rsym.w_m --> roughly(rsym.w_f)
@fact rsym.u_m --> roughly(rsym.u_f)
@fact rsym.α --> roughly(rsym.α')

# asymmetric case: needs σ >~ 10 to converge
rasym = closed_uniform(50, 50, 100, h, σ)

# check convergence
rmsse, rfsse = sse_resid(rasym)
rmvf, rfvf = vf_resid(rasym)

# valid solution
@fact rmsse --> roughly(zeros(rmsse), atol = 1e-7)
@fact rfsse --> roughly(zeros(rfsse), atol = 1e-7)
# convergence seems to stall around 5e-6...
@fact maximum(abs.(rmvf)) --> roughly(0.0, atol = 1e-5)
@fact maximum(abs.(rfvf)) --> roughly(0.0, atol = 1e-5)
@fact rasym.α --> 1.0 .- G.(-match_surplus(rasym))

# sex ratio effects on singles
@fact (rsym.u_f .≤ rasym.u_f) --> all

# inflow model: symmetric
inflow_symm = inflow_uniform(50, 100, 100, h, σ)
imsse, ifsse = sse_resid(inflow_symm)
@fact inflow_symm.w_m --> roughly(inflow_symm.w_f)
@fact inflow_symm.u_m --> roughly(inflow_symm.u_f)
@fact inflow_symm.α --> roughly(inflow_symm.α')
@fact imsse --> roughly(zeros(imsse), atol = 1e-7)
@fact ifsse --> roughly(zeros(ifsse), atol = 1e-7)

# inflow model: asymmetric
inflow_asymm = inflow_uniform(50, 50, 100, h, σ)
rimsse, rifsse = sse_resid(inflow_asymm)
@fact rimsse --> roughly(zeros(rimsse), atol = 1e-7)
@fact rifsse --> roughly(zeros(rifsse), atol = 1e-7)
