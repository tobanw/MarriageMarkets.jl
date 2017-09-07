# uses variables from setup_tests.jl

STDNORMAL = Normal()
G(x::Real) = cdf(Normal(0, σ), x)
μ(a::Real) = σ * (pdf(STDNORMAL, quantile(STDNORMAL, 1-a)) - a * quantile(STDNORMAL, 1-a))


### Verifications: endogenous divorce model ###

"""
Element-by-element calculation of steady state flow equations
```julia
mres = ℓ_m - um .* (1 .+ λ .* ((α ./ (δ*(1-α) .+ ψ_m .+ ψ_f')) * uf))
fres = ℓ_f - uf .* (1 .+ λ .* ((α ./ (δ*(1-α) .+ ψ_m .+ ψ_f'))' * um))
```
"""
function sse_stoch(M::SearchMatch)

	mres = similar(M.ℓ_m)
	fres = similar(M.ℓ_f)
	
	for i in 1:length(M.ℓ_m)
		mres[i] = M.ℓ_m[i] - M.u_m[i] * (1 + M.λ * dot(M.α[i,:] ./
		          (M.δ * (1 - M.α[i,:]) + M.ψ_m[i] + M.ψ_f), M.u_f))
	end
	for j in 1:length(M.ℓ_f)
		fres[j] = M.ℓ_f[j] - M.u_f[j] * (1 + M.λ * dot(M.α[:,j] ./
		          (M.δ * (1 - M.α[:,j]) + M.ψ_f[j] + M.ψ_m), M.u_m))
	end
	return mres, fres
end

"""
Element-by-element calculation of value function equations
```julia

v_m = 0.5*λ * (μα * u_f)
v_f = 0.5*λ * (μα' * u_m)
```
where `μα = μ.(α) ./ (r + δ + ψ_m .+ ψ_f')`
"""
function vf_stoch(M::SearchMatch)
	mres = similar(M.v_m)
	fres = similar(M.v_f)

	for i in 1:length(M.v_m)
		mres[i] = 2*M.v_m[i] - M.λ * dot(μ.(M.α[i,:]) ./
		                             (M.r + M.δ + M.ψ_m[i] + M.ψ_f), M.u_f)
	end
	for j in 1:length(M.v_f)
		fres[j] = 2*M.v_f[j] - M.λ * dot(μ.(M.α[:,j]) ./
		                             (M.r + M.δ + M.ψ_f[j] + M.ψ_m), M.u_m)
	end

	return mres, fres
end

function surplus_stoch(M::SearchMatch)
	s = similar(M.h)

	for i in 1:length(M.v_m), j in 1:length(M.v_f)
		s[i,j] = M.h[i,j] - M.v_m[i] - M.v_f[j] + M.δ * μ(M.α[i,j]) /
		                               (M.r + M.δ + M.ψ_m[i] + M.ψ_f[j])
	end

	return s
end


### Tests: endogenous divorce model ###

# symmetric case
rsym = fetch(rsym_job) # get result from worker process

@fact rsym.v_m --> roughly(rsym.v_f) "expected symmetry of values"
@fact rsym.u_m --> roughly(rsym.u_f) "expected symmetry of singles"
@fact rsym.α --> roughly(rsym.α') "expected symmetry of matching"


# asymmetric case
rasym = fetch(rasym_job) # get result from worker process

# check convergence
rmsse, rfsse = sse_stoch(rasym)
rmvf, rfvf = vf_stoch(rasym)

# valid solution
@fact rmsse --> roughly(zeros(rmsse), atol=1e-7) "market equilibrium: single men did not converge"
@fact rfsse --> roughly(zeros(rfsse), atol=1e-7) "market equilibrium: single women did not converge"
# convergence seems to stall around 2e-5...
@fact maximum(abs.(rmvf)) --> roughly(0, atol=3e-5) "matching equilibrium: man values did not converge"
@fact maximum(abs.(rfvf)) --> roughly(0, atol=1e-7) "matching equilibrium: woman values did not converge"
@fact rasym.α --> roughly(1 - G.(-surplus_stoch(rasym)), atol=1e-5) "matching equilibrium: match probabilities α did not converge"

# sex ratio effects on singles
@fact (rsym.u_f .≤ rasym.u_f) --> all "expected more singlehood with fewer available men"
