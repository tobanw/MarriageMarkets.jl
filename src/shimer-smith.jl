using NLsolve
using QuantEcon: compute_fixed_point

"""
Construct a Shimer & Smith (2000) marriage market model and solve for the equilibrium.

The equilibrium is the solution to a three-part fixed-point mapping.
Without randomness, the fixed point iteration will not converge with any sex asymmetry.
"""
type ShimerSmith

	### Parameters ###

	"arrival rate of meetings"
	ρ::Float64
	"arrival rate of separations"
	δ::Float64
	"discount rate"
	r::Float64


	### Exogenous objects ###

	"male population distribution: not normalized"
	ℓ_m::Vector
	"female population distribution: not normalized"
	ℓ_f::Vector

	"production function as array"
	h::Array

	"CDF of match-specific additive shock"
	G::Function


	### Endogenous equilibrium objects ###

	"mass of single males"
	u_m::Vector
	"mass of single females"
	u_f::Vector

	"match function as array"
	α::Array

	"male singlehood (average) value function as vector"
	w_m::Vector
	"female singlehood (average) value function as vector"
	w_f::Vector

	"marital surplus function as array"
	S::Array


	### Inner Constructor ###

	"""
		ShimerSmith(ρ, δ, r, ℓ_m, ℓ_f, h, G)

	Solve the model with match-specific shocks ``z ~ G(z)``.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64, ℓ_m::Vector, ℓ_f::Vector,
					  h::Array, G::Function)

		# argument validition
		if minimum(ℓ_m) < 0.0 || minimum(ℓ_f) < 0.0
			error("Population masses must be non-negative.")
		elseif length(ℓ_m) != size(h)[1] || length(ℓ_f) != size(h)[2]
			error("Number of types inconsistent with production array.")
		elseif G(Inf) ≠ 1.0 || G(-Inf) ≠ 0.0 ||
			any(G.(linspace(-10.5, 9.5, 100)) .> G.(linspace(-9.5, 10.5, 100)))
			error("G not a valid CDF.")
		end


		# Initialize guesses for w: overwritten and reused in the inner fixed point iteration
		w_m = 0.5 * h[:,1]
		w_f = 0.5 * h[1,:]

		# Initialize matching array: overwritten and reused in the outer fixed point iteration
		α = ones(Float64, length(ℓ_m), length(ℓ_f))

		"""
		Matching equilibrium fixed point operator ``T_α(α)``.
		"""
		function fp_matching_eqm(A::Array, u_m::Vector, u_f::Vector)
			# overwrite w_m, w_f to reuse as initial guess for nlsolve
			w_m[:], w_f[:] = update_values(ρ, δ, r, A, h, u_m, u_f, w_m, w_f)
			return update_match(G, h, w_m, w_f)
		end

		"""
		Market equilibrium fixed point operator ``T_u(u_m, u_f)``.
		"""
		function fp_market_eqm(u::Vector; inner_tol=1e-4)

			um, uf = sex_split(u, length(ℓ_m))

			# nested fixed point: overwrite α to be reused as initial guess for next call
			α[:] = compute_fixed_point(x->fp_matching_eqm(x, um, uf), α,
								 err_tol=inner_tol, verbose=1)

			# steady state distributions
			um_new, uf_new = update_singles(ρ, δ, α, ℓ_m, ℓ_f, um, uf)

			# truncate if `u` strays out of bounds
			if minimum([um_new; uf_new]) < 0.0
				warn("u negative: truncating...")
			elseif minimum([ℓ_m .- um_new; ℓ_f .- uf_new]) < 0.0
				warn("u > ℓ: truncating...")
			end
			um_new[:] = clamp.(um_new, 0.0, ℓ_m)
			uf_new[:] = clamp.(uf_new, 0.0, ℓ_f)

			return [um_new; uf_new]
		end

		# fast rough compututation of equilibrium by fixed point iteration
		u_fp0 = compute_fixed_point(fp_market_eqm, 0.1*[ℓ_m; ℓ_f],
							  print_skip=1, verbose=2) # initial guess u = 0.1*ℓ

		um_0, uf_0 = sex_split(u_fp0, length(ℓ_m))

		# touch up with high precision fixed point solution
		u_fp = compute_fixed_point(x->fp_market_eqm(x, inner_tol=1e-10), [um_0; uf_0],
							 err_tol=1e-10, print_skip=5, verbose=1)
		
		u_m, u_f = sex_split(u_fp, length(ℓ_m))

		S = match_surplus(h, w_m, w_f)

		# construct instance
		new(ρ, δ, r, ℓ_m, ℓ_f, h, G, u_m, u_f, α, w_m, w_f, S)

	end # constructor


	### Outer Constructors ###

	"""
		ShimerSmith(ρ, δ, r, Θ_m, Θ_f, ℓ_m, ℓ_f, g, G)

	Solve the model with match-specific shocks ``z ~ G(z)`` and production function ``g(x,y)``.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64,
					  Θ_m::Vector, Θ_f::Vector, ℓ_m::Vector, ℓ_f::Vector,
					  g::Function, G::Function)
		h = prod_array(Θ_m, Θ_f, g)
		return ShimerSmith(ρ, δ, r, ℓ_m, ℓ_f, h, G)
	end

	"""
		ShimerSmith(ρ, δ, r, Θ_m, Θ_f, ℓ_m, ℓ_f, g)

	Solve the model without randomness and with production function ``g(x,y)``.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64,
					  Θ_m::Vector, Θ_f::Vector, ℓ_m::Vector, ℓ_f::Vector, g::Function)
		# degenerate distribution: point mass of one at x=0, but no marriage if indifferent
		G(x) = Float64(x > 0.0) # bool as float
		return ShimerSmith(ρ, δ, r, Θ_m, Θ_f, ℓ_m, ℓ_f, g, G)
	end


	### Equilibrium computation ###

	"""
		update_singles(ρ, δ, α, ℓ_m, ℓ_f, u_m, u_f)

	Compute the implied singles distributions given a matching function from the
	steady-state equilibrium conditions.

	For males, this equation is:
	``∀x, δ(ℓ(x) - u(x)) = ρ u(x) ∫ α(x,y) u(y) dy``.
	Thus, this function solves a non-linear system of equations for `u_m` and `u_f`.

	The constraints ``0 ≤ u ≤ ℓ`` are not enforced here, but the outputs of this function
	are truncated in the fixed point iteration loop.
	"""
	function update_singles(ρ::Float64, δ::Float64, α::Array,
						 ℓ_m::Vector, ℓ_f::Vector, u_m::Vector, u_f::Vector)

		"""
		Vector function representing the system of equations to be solved.
		"""
		function vecsys!(μ::Vector, res::Vector) # vec'd
			μm, μf = sex_split(μ, length(ℓ_m))

			# compute residuals of non-linear system
			mres = δ * ℓ_m - μm .* (δ + ρ * (α * μf))
			fres = δ * ℓ_f - μf .* (δ + ρ * (α' * μm))

			res[:] = [mres; fres] # concatenate into stacked vector
		end # vecsys!

		# initial guess of shares of singles: stacked vector of previous values
		guess = [u_m; u_f]

		# NLsolve
		result = nlsolve(vecsys!, guess)

		um_new, uf_new = sex_split(result.zero, length(ℓ_m))

		return um_new, uf_new

	end # update_singles

	"""
		update_values(ρ, δ, r, α, h, u_m, u_f, w_m, w_f)

	Compute the implied singlehood value functions given a matching function
	and singles distributions.

	For males, the value function equation is:
	``∀x, w(x) = \frac{ρ}{2(r+δ)} ∫ (h(x,y) - w(x) - w(y)) α(x,y) u(y) dy``.
	Thus, this function solves a non-linear system of equations for `w_m` and `w_f`.
	"""
	function update_values(ρ::Float64, δ::Float64, r::Float64, α::Array, h::Array,
						 u_m::Vector, u_f::Vector, w_m::Array, w_f::Array)
		θ = ρ / (2*(r+δ)) # precompute constant

		"""
		Vector function representing the system of equations to be solved.
		"""
		function vecsys!(ω::Vector, res::Vector) # vec'd
			ωm, ωf = sex_split(ω, length(w_m))

			αS = α .* match_surplus(h, ωm, ωf)

			# compute residuals of non-linear system
			mres = ωm - θ * (αS * u_f)
			fres = ωf - θ * (αS' * u_m)

			res[:] = [mres; fres] # concatenate into stacked vector
		end # vecsys!

		# initial guess of value functions: stacked vector of previous values
		guess = [w_m; w_f]

		# NLsolve
		result = nlsolve(vecsys!, guess)

		wm_new, wf_new = sex_split(result.zero, length(w_m))

		return wm_new, wf_new

	end # update_values

	"""
	Calculate matching function ``α(x,y)`` from ``S ≥ 0`` condition.

	When `G` is degenerate, this yields the non-random case.
	"""
	function update_match(G::Function, h::Array, w_m::Array, w_f::Array)
		return 1.0 .- G.(-match_surplus(h, w_m, w_f))
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

	"Construct match surplus array from value functions."
	function match_surplus(h::Array, w_m::Vector, w_f::Vector)
		S = similar(h)

		for i in 1:length(w_m), j in 1:length(w_f)
			S[i,j] = h[i,j] - w_m[i] - w_f[j]
		end

		return S
	end # match_surplus

	"Split vector `v` into male/female pieces, where `idx` is number of male types."
	function sex_split(v::Vector, idx::Int)
		vm = v[1:idx]
		vf = v[idx+1:end]
		return vm, vf
	end

end # type
