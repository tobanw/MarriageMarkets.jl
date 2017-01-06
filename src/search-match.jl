using NLsolve
using QuantEcon: compute_fixed_point
using Distributions

const STDNORMAL = Normal()

"""
	SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h)

Construct a Shimer & Smith (2000) marriage market model and solve for the equilibrium.

The equilibrium is the solution to a three-part fixed-point mapping.

Model selection depends which arguments are provided:
* Match-specific additive shocks ``z ~ N(0, σ)``
	* Note: some randomness is required for the fixed point iteration to converge with sex asymmetry
* Closed system or inflow/outflow:
	* `ℓ_m, ℓ_f` exogenous: population cycles between singlehood and marriage, no birth/death
	* Death rates `ψ_m, ψ_f`, inflows `γ_m, γ_f`: population distributions `ℓ_m, ℓ_f` endogenous
"""
immutable SearchMatch # object fields cannot be modified

	### Parameters ###

	"arrival rate of meetings"
	ρ::Real
	"arrival rate of separation shocks"
	δ::Real
	"discount rate"
	r::Real
	"standard deviation of normally distributed match-specific additive shock"
	σ::Real

	### Exogenous objects ###

	"male inflows by type"
	γ_m::Vector
	"female inflows by type"
	γ_f::Vector

	"arrival rates of male death by type"
	ψ_m::Vector
	"arrival rates of female death by type"
	ψ_f::Vector

	"production function as array"
	h::Array


	### Endogenous equilibrium objects ###

	"male population distribution: not normalized"
	ℓ_m::Vector
	"female population distribution: not normalized"
	ℓ_f::Vector

	"mass of single males"
	u_m::Vector
	"mass of single females"
	u_f::Vector

	"male singlehood (average) value function as vector"
	w_m::Vector
	"female singlehood (average) value function as vector"
	w_f::Vector

	"match function as array"
	α::Array

	"marital surplus function as array"
	S::Array


	### Inner Constructor ###

	"""
	Inner constructor that accomodates several model variants, depending on the arguments provided.
	This is not meant to be called directly -- instead, outer constructors should call this
	constructor with a full set of arguments, using zero (or identity) values for unwanted components.
	"""
	#TODO: change G to σ
	function SearchMatch(ρ::Real, δ::Real, r::Real, σ::Real,
                      γ_m::Vector, γ_f::Vector, ψ_m::Vector, ψ_f::Vector,
                      ℓ_m::Vector, ℓ_f::Vector, h::Array)

		### Model Selection ###

		const n_m = length(γ_m)
		const n_f = length(γ_f)

		if sum(ψ_m) > 0 && sum(ψ_f) > 0 # inflow/outflow model: if death rates provided
			const INFLOW = true
		else
			const INFLOW = false
		end

		if σ == 0
			const STOCH = false
		else
			const STOCH = true
		end

		const AGING = false # TODO


		### Argument Validation ###

		if INFLOW # birth/death model
			if any(ℓ_m .> 0) || any(ℓ_f .> 0)
				error("Death rate ψ provided: population distributions ℓ_m, ℓ_f are endogenous.")
			elseif sum(γ_m) ≤ 0 || sum(γ_f) ≤ 0
				error("Total population inflow must be positive.")
			elseif any(ψ_m .< 0) || any(ψ_f .< 0)
				error("Death rates must be non-negative.")
			end
		else # closed system: population exogenous
			if sum(ℓ_m) ≤ 0 || sum(ℓ_f) ≤ 0
				error("No death: population must be positive.")
			elseif any(ℓ_m .< 0) || any(ℓ_f .< 0)
				error("Population masses must be non-negative.")
			elseif length(ℓ_m) != size(h)[1] || length(ℓ_f) != size(h)[2]
				error("Number of types inconsistent with production array.")
			end
		end

		if STOCH && σ < 0
			error("σ must be non-negative.")
		end

		# more argument validation
		if any([ρ, δ, r] .≤ 0)
			error("Parameters ρ, δ, r must be positive.")
		elseif size(h)[1] != n_m || size(h)[2] != n_f
			error("Number of types inconsistent with production array.")
		elseif length(ℓ_m) ≠ n_m || length(ℓ_f) ≠ n_f
			error("Inconsistent number of types.")
		end


		### Compute Constants and Values ###

		if INFLOW
			# population size given directly by inflow and outflow rates
			ℓ_m = γ_m ./ ψ_m
			ℓ_f = γ_f ./ ψ_f
		end

		# cdf of match-specific marital productivity shocks
		G(x::Real) = STOCH ? cdf(Normal(0, σ), x) : Float64(x ≥ 0) # bool as float


		### Steady-State Equilibrium Conditions ###
		"""
		Update population distributions: basic Shimer-Smith model.

		Compute the implied singles distributions given a matching function.
		For males, the steady state condition equating flows into and out of marriage is:
		```math
		∀x, (δ + ψ(x))(ℓ(x) - u(x)) = ρ u(x) ∫ α(x,y) u(y) dy
		```
		Thus, this function solves a non-linear system of equations for `u_m` and `u_f`.
		The constraints ``0 ≤ u ≤ ℓ`` are not enforced here, but the outputs of this function
		are truncated in the fixed point iteration loop.
		"""
		#TODO: update with endogenous divorce
		function steadystate_base!(u::Vector, res::Vector) # stacked vector
			# uses the overwritable α in the outer scope

			um, uf = sex_split(u, n_m)

			# compute residuals of non-linear system
			mres = (δ .+ ψ_m) .* ℓ_m - um .* ((δ .+ ψ_m) .+ ρ .* (α * uf))
			fres = (δ .+ ψ_f) .* ℓ_f - uf .* ((δ .+ ψ_f) .+ ρ .* (α' * um))

			res[:] = [mres; fres] # concatenate into stacked vector
		end # steadystate_base!

		"Update population distributions: aging model."
		function steadystate_aging!(u::Vector, res::Vector) # stacked vector
			# uses the overwritable α in the outer scope

			um, uf = sex_split(u, n_m)

			#TODO: add actual conditions
			# compute residuals of non-linear system
			mres = (δ .+ ψ_m) .* ℓ_m - um .* ((δ .+ ψ_m) .+ ρ .* (α * uf))
			fres = (δ .+ ψ_f) .* ℓ_f - uf .* ((δ .+ ψ_f) .+ ρ .* (α' * um))

			res[:] = [mres; fres] # concatenate into stacked vector
		end # steadystate_inflow!

		"""
		Update singlehood value functions: basic Shimer-Smith model.

		Compute the implied average singlehood value functions given a matching function
		and singles distributions.
		For males, the average value function equation is:
		```math
		∀x, 2w(x) = ρ ∫ α(x,y) S(x,y) u(y) dy,\\
		S(x,y) = \frac{h(x,y) - w(x) - w(y)}{r+δ+ψ_m(x)+ψ_f(y)}
		```
		Thus, this function solves a non-linear system of equations for `w_m` and `w_f`,
		the average value functions, which are defined in terms of their present-value counterparts as
		```math
		w_m(x) = (r+ψ_m(x))W_m(x)
		```
		"""
		#TODO: update with match shocks and endogenous divorce
		# if STOCH: EV of shock is mean of truncated normal, exp_trunc(a,σ)
		# else: use s without EV term
		function valuefunc_base!(ω::Vector, res::Vector, u_m::Vector, u_f::Vector) # stacked vector
			ωm, ωf = sex_split(ω, n_m)

			αS = α .* match_surplus(ωm, ωf)

			# compute residuals of non-linear system
			mres = 2*ωm - ρ * (αS * u_f)
			fres = 2*ωf - ρ * (αS' * u_m)

			res[:] = [mres; fres] # concatenate into stacked vector
		end # vecsys!

		"""
		Calculate matching function ``α(x,y)`` from ``S ≥ 0`` condition.

		When `G` is degenerate, this yields the non-random case.
		"""
		function update_match(w_m::Array, w_f::Array)
			return 1.0 .- G.(-match_surplus(w_m, w_f))
		end

		"Construct match surplus array S from value functions."
		# TODO: implement surplus for divorce
		function match_surplus(w_m::Vector, w_f::Vector)
			S = similar(h)

			for i in 1:length(w_m), j in 1:length(w_f)
				S[i,j] = (h[i,j] - w_m[i] - w_f[j]) / (r + δ + ψ_m[i] + ψ_f[j])
			end

			return S
		end # match_surplus


		### Equilibrium Solver ###

		# Initialize guesses for w: overwritten and reused in the inner fixed point iteration
		w_m = 0.5 * h[:,1]
		w_f = 0.5 * h[1,:]

		# Initialize matching array: overwritten and reused in the outer fixed point iteration
		α = ones(Float64, n_m, n_f)

		"""
		Matching equilibrium fixed point operator ``T_α(α)``.
		"""
		function fp_matching_eqm(A::Array, u_m::Vector, u_f::Vector)
			# overwrite w_m, w_f to reuse as initial guess for nlsolve
			w_m[:], w_f[:] = sex_solve((x,res)->valuefunc_base!(x, res, u_m, u_f), w_m, w_f)
			return update_match(w_m, w_f)
		end

		"""
		Market equilibrium fixed point operator ``T_u(u_m, u_f)``.
		"""
		function fp_market_eqm(u::Vector; inner_tol=1e-4)

			um, uf = sex_split(u, n_m)

			# nested fixed point: overwrite α to be reused as initial guess for next call
			α[:] = compute_fixed_point(x->fp_matching_eqm(x, um, uf), α,
								 err_tol=inner_tol, verbose=1)

			# steady state distributions
			if AGING
				warn("Aging model not yet implemented.")
				um_new, uf_new = sex_solve(steadystate_aging!, um, uf)
			else
				um_new, uf_new = sex_solve(steadystate_base!, um, uf)
			end

			# truncate if `u` strays out of bounds
			if minimum([um_new; uf_new]) < 0
				warn("u negative: truncating...")
			elseif minimum([ℓ_m .- um_new; ℓ_f .- uf_new]) < 0
				warn("u > ℓ: truncating...")
			end
			um_new[:] = clamp.(um_new, 0, ℓ_m)
			uf_new[:] = clamp.(uf_new, 0, ℓ_f)

			return [um_new; uf_new]
		end

		# fast rough compututation of equilibrium by fixed point iteration
		u_fp0 = compute_fixed_point(fp_market_eqm, 0.1*[ℓ_m; ℓ_f],
                              print_skip=10, verbose=2) # initial guess u = 0.1*ℓ

		um_0, uf_0 = sex_split(u_fp0, n_m)

		# touch up with high precision fixed point solution
		u_fp = compute_fixed_point(x->fp_market_eqm(x, inner_tol=1e-5), [um_0; uf_0],
							 err_tol=1e-8, verbose=1)
		
		u_m, u_f = sex_split(u_fp, n_m)

		S = match_surplus(w_m, w_f)

		# construct instance
		new(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, h, ℓ_m, ℓ_f, u_m, u_f, w_m, w_f, α, S)

	end # constructor

end # type


### Outer Constructors ###

# Recall inner constructor: SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h)

"Closed-system model with match-specific gaussian shocks and production function ``g(x,y)``."
function SearchClosed(ρ::Real, δ::Real, r::Real, σ::Real,
                      Θ_m::Vector, Θ_f::Vector, ℓ_m::Vector, ℓ_f::Vector,
                      g::Function)
	# irrelevant arguments to pass as zeros
	ψ_m = zeros(ℓ_m)
	ψ_f = zeros(ℓ_f)
	γ_m = zeros(ℓ_m)
	γ_f = zeros(ℓ_f)

	h = prod_array(Θ_m, Θ_f, g)
	return SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h)
end

"Inflow model with match-specific gaussian shocks and production function ``g(x,y)``."
function SearchInflow(ρ::Real, δ::Real, r::Real, σ::Real,
                      Θ_m::Vector, Θ_f::Vector, γ_m::Vector, γ_f::Vector,
                      ψ_m::Vector, ψ_f::Vector, g::Function)
	# irrelevant arguments to pass as zeros
	ℓ_m = zeros(γ_m)
	ℓ_f = zeros(γ_f)

	h = prod_array(Θ_m, Θ_f, g)
	return SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h)
end


### Helper functions ###

"Construct production array from function."
function prod_array(mtypes::Vector, ftypes::Vector, prodfn::Function)
	h = Array{Float64}(length(mtypes), length(ftypes))

	for (i,x) in enumerate(mtypes), (j,y) in enumerate(ftypes)
		h[i,j] = prodfn(x, y)
	end

	return h

end # prod_array

"Wrapper for nlsolve that handles the concatenation and splitting of sex vectors."
function sex_solve(eqnsys!, v_m, v_f)
	# initial guess: stacked vector of previous values
	guess = [v_m; v_f]

	# NLsolve
	result = nlsolve(eqnsys!, guess)

	vm_new, vf_new = sex_split(result.zero, length(v_m))

	return vm_new, vf_new
end

"Split vector `v` into male/female pieces, where `idx` is number of male types."
function sex_split(v::Vector, idx::Int)
	vm = v[1:idx]
	vf = v[idx+1:end]
	return vm, vf
end

"Upper tail truncated normal: E[z|z>a] = σ ϕ(a/σ) / (1 - Φ(a/σ))."
function exp_trunc(a::Real, σ::Real)
	return σ * pdf(STDNORMAL, a/σ) / (1 - cdf(STDNORMAL, a/σ))
end
