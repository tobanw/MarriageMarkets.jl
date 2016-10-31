ρ = 5.0
r = 0.05
δ = 0.05

function SS_uniform(ntypes, mmass, fmass, prod)
	
	# types
	Θ = Vector(linspace(0.0, 1.0, ntypes))

	# uniform population distributions
	lm = (mmass / ntypes) * ones(Float64, ntypes)
	lf = (fmass / ntypes) * ones(Float64, ntypes)

	return ShimerSmith(ρ, δ, r, Θ, Θ, lm, lf, prod)
end

h(x::Real, y::Real) = x*y

function sse_resid(M::ShimerSmith)
	#=
	element-by-element calculation of
	mres = M.δ * M.ℓ_m - M.u_m .* (M.δ + M.ρ * (M.α * M.u_f))
	fres = M.δ * M.ℓ_f - M.u_f .* (M.δ + M.ρ * (M.α' * M.u_m))
	=#

	mres = similar(M.ℓ_m)
	fres = similar(M.ℓ_f)
	
	for i in 1:length(M.ℓ_m)
		mres[i] = M.δ * M.ℓ_m[i] - M.u_m[i] * (M.δ + M.ρ * dot(M.α[i,:], M.u_f))
	end
	for j in 1:length(M.ℓ_f)
		fres[j] = M.δ * M.ℓ_f[j] - M.u_f[j] * (M.δ + M.ρ * dot(M.α[:,j], M.u_m))
	end
	return mres, fres
end

function vf_resid(M::ShimerSmith)
	#=
	element-by-element calculation of
	mres = M.w_m - θ * (αS * M.u_f)
	fres = M.w_f - θ * (αS' * M.u_m),
	where:
	αS = M.α .* M.S
	=#

	θ = M.ρ / (2*(M.r+M.δ)) # precompute constant

	mres = similar(M.w_m)
	fres = similar(M.w_f)

	for i in 1:length(M.w_m)
		mres[i] = M.w_m[i] - θ * dot(M.α[i,:] .* M.S[i,:], M.u_f)
	end
	for j in 1:length(M.w_f)
		fres[j] = M.w_f[j] - θ * dot(M.α[:,j] .* M.S[:,j], M.u_m)
	end
	return mres, fres
end

function match_surplus(h::Array, w_m::Vector, w_f::Vector)
	S = similar(h)

	for i in 1:length(w_m), j in 1:length(w_f)
		S[i,j] = h[i,j] - w_m[i] - w_f[j]
	end

	return S
end

# symmetric case
symm = SS_uniform(50, 100, 100, h)

# check convergence
msse, fsse = sse_resid(symm)
mvf, fvf = vf_resid(symm)

# valid solution
@fact msse --> roughly(zeros(msse), atol = 1e-8)
@fact fsse --> roughly(zeros(fsse), atol = 1e-8)
@fact mvf --> roughly(zeros(mvf), atol = 1e-8)
@fact fvf --> roughly(zeros(fvf), atol = 1e-8)
@fact symm.α --> 1*(match_surplus(symm.h, symm.w_m, symm.w_f) .≥ 0.0)

# symmetry
@fact symm.w_m --> roughly(symm.w_f)
@fact symm.u_m --> roughly(symm.u_f)
@fact symm.α --> roughly(symm.α')

#= FIXME: asymmetric case doesn't converge!
asymm = SS_uniform(50, 50, 100, h)

# sex ratio effects on value functions
@fact (symm.w_f .≥ asymm.w_f) --> all
=#
