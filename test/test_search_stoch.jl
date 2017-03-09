### Helper functions ###

using Distributions

ρ = 5.0
r = 0.05
δ = 0.05
σ = 10

const STDNORMAL = Normal()
μ(a::Real) = σ * (pdf(STDNORMAL, quantile(STDNORMAL, 1-a)) - a * quantile(STDNORMAL, 1-a))

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


### Verifications: endogenous divorce model ###

"""
Element-by-element calculation of steady state flow equations
```julia
mres = ℓ_m - um .* (1 .+ ρ .* ((α ./ (δ*(1-α) .+ ψ_m .+ ψ_f')) * uf))
fres = ℓ_f - uf .* (1 .+ ρ .* ((α ./ (δ*(1-α) .+ ψ_m .+ ψ_f'))' * um))
```
"""
function sse_stoch(M::SearchMatch)

	mres = similar(M.ℓ_m)
	fres = similar(M.ℓ_f)
	
	for i in 1:length(M.ℓ_m)
		mres[i] = M.ℓ_m[i] - M.u_m[i] * (1 + M.ρ * dot(M.α[i,:] ./
		          (M.δ * (1 - M.α[i,:]) + M.ψ_m[i] .+ M.ψ_f), M.u_f))
	end
	for j in 1:length(M.ℓ_f)
		fres[j] = M.ℓ_f[j] - M.u_f[j] * (1 + M.ρ * dot(M.α[:,j] ./
		          (M.δ * (1 - M.α[:,j]) + M.ψ_f[j] .+ M.ψ_m), M.u_m))
	end
	return mres, fres
end

"""
Element-by-element calculation of value function equations
```julia

v_m = 0.5*ρ * (μα * u_f)
v_f = 0.5*ρ * (μα' * u_m)
```
where `μα = μ.(α) ./ (r + δ + ψ_m .+ ψ_f')`
"""
function vf_stoch(M::SearchMatch)
	mres = similar(M.v_m)
	fres = similar(M.v_f)

	for i in 1:length(M.v_m)
		mres[i] = 2*M.v_m[i] - M.ρ * dot(μ.(M.α[i,:]) ./
		                             (M.r + M.δ + M.ψ_m[i] .+ M.ψ_f), M.u_f)
	end
	for j in 1:length(M.v_f)
		fres[j] = 2*M.v_f[j] - M.ρ * dot(μ.(M.α[:,j]) ./
		                             (M.r + M.δ + M.ψ_f[j] .+ M.ψ_m), M.u_m)
	end

	return mres, fres
end

function surplus_stoch(M::SearchMatch)
	s = similar(M.h)

	for i in 1:length(M.v_m), j in 1:length(M.v_f)
		s[i,j] = M.h[i,j] - M.v_m[i] - M.v_f[j] + M.δ * μ(M.α[i,j]) /
		                               (M.r + M.δ + M.ψ_m[i] .+ M.ψ_f[j])
	end

	return s
end


### Tests: endogenous divorce model ###

h(x::Real, y::Real) = x*y
G(x::Real) = cdf(Normal(0, σ), x)

# symmetric case
rsym = closed_uniform(50, 100, 100, h, σ)

@fact rsym.v_m --> roughly(rsym.v_f)
@fact rsym.u_m --> roughly(rsym.u_f)
@fact rsym.α --> roughly(rsym.α')

# asymmetric case: needs σ >~ 10 to converge
rasym = closed_uniform(50, 50, 100, h, σ)

# check convergence
rmsse, rfsse = sse_stoch(rasym)
rmvf, rfvf = vf_stoch(rasym)

# valid solution
@fact rmsse --> roughly(zeros(rmsse), atol = 1e-7)
@fact rfsse --> roughly(zeros(rfsse), atol = 1e-7)
# convergence seems to stall around 2e-5...
@fact maximum(abs.(rmvf)) --> roughly(0.0, atol = 3e-5)
@fact maximum(abs.(rfvf)) --> roughly(0.0, atol = 1e-7)
@fact rasym.α --> roughly(1.0 .- G.(-surplus_stoch(rasym)), atol=1e-5)

# sex ratio effects on singles
@fact (rsym.u_f .≤ rasym.u_f) --> all

# inflow model: symmetric
inflow_symm = inflow_uniform(50, 100, 100, h, σ)
imsse, ifsse = sse_stoch(inflow_symm)
@fact inflow_symm.v_m --> roughly(inflow_symm.v_f)
@fact inflow_symm.u_m --> roughly(inflow_symm.u_f)
@fact inflow_symm.α --> roughly(inflow_symm.α')
@fact imsse --> roughly(zeros(imsse), atol = 1e-7)
@fact ifsse --> roughly(zeros(ifsse), atol = 1e-7)

# inflow model: asymmetric
inflow_asymm = inflow_uniform(50, 50, 100, h, σ)
rimsse, rifsse = sse_stoch(inflow_asymm)
@fact rimsse --> roughly(zeros(rimsse), atol = 1e-7)
@fact rifsse --> roughly(zeros(rifsse), atol = 1e-7)
