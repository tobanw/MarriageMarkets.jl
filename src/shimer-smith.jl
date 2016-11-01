using NLsolve
using QuantEcon: compute_fixed_point

"""
Construct a Shimer & Smith (2000) marriage market model and solve for the equilibrium.

The equilibrium is the solution to a three-part fixed-point mapping.
Without randomness, the fixed point iteration will not converge with any sex asymmetry.
"""
type ShimerSmith

	### Parameters ###

	ρ::Float64 # arrival rate of meetings
	δ::Float64 # arrival rate of separations
	r::Float64 # discount rate


	### Exogenous objects ###

	# population distributions: not normalized
	ℓ_m::Vector
	ℓ_f::Vector

	# production function as array
	h::Array

	# match-specific additive shock CDF
	G::Function


	### Endogenous equilibrium objects ###

	# masses of singles
	u_m::Vector
	u_f::Vector

	# match function as array
	α::Array

	# (average) value functions of singles as vectors
	w_m::Vector
	w_f::Vector

	# marital surplus function as array
	S::Array


	### Inner Constructor ###

	"""
		ShimerSmith(ρ, δ, r, ℓ_m, ℓ_f, h, G)

	Solve the model with match-specific shocks ``z ~ G(z)``.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64, ℓ_m::Vector, ℓ_f::Vector,
					  h::Array, G::Function)

		# CHECK: masses of m/f must be valid distro
		if minimum(ℓ_m) < 0.0 || minimum(ℓ_f) < 0.0
			error("Invalid type distribution.")
		end

		# CHECK: masses of m/f match the number of types
		if length(ℓ_m) != size(h)[1] || length(ℓ_f) != size(h)[2]
			error("Number of types inconsistent with production array.")
		end

		# Initialize objects: guesses
		u_m = 1.0 * ℓ_m
		u_f = 1.0 * ℓ_f
		w_m = 0.5 * h[:,1]
		w_f = 0.5 * h[1,:]
		a = ones(Float64, length(ℓ_m), length(ℓ_f))

		"""
		Fixed point operator ``T(w)``.
		"""
		function fp_operator(A::Array)
			# compose mappings
			update_singles!(ρ, δ, A, ℓ_m, ℓ_f, u_m, u_f) # update u_m, u_f in place

			# truncate if `u` strays out of bounds
			u_m[:] = clamp.(u_m, 0.0, ℓ_m)
			u_f[:] = clamp.(u_f, 0.0, ℓ_f)

			update_values!(ρ, δ, r, A, h, u_m, u_f, w_m, w_f) # update w_m, w_f in place

			return update_match(G, h, w_m, w_f)
		end

		# compute the equilibrium by fixed point iteration
		α = compute_fixed_point(fp_operator, a, err_tol=1e-10, print_skip=1, verbose=2)

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
		# degenerate distribution: point mass of one at x=0
		G(x) = 1.0 * (x ≥ 0.0)
		return ShimerSmith(ρ, δ, r, Θ_m, Θ_f, ℓ_m, ℓ_f, g, G)
	end


	### Equilibrium computation ###

	"""
		update_singles!(ρ, δ, α, ℓ_m, ℓ_f, u_m, u_f)

	Compute the implied singles distributions given a matching function from the
	steady-state equilibrium conditions.

	For males, this equation is:
	``∀x, δ(ℓ(x) - u(x)) = ρ u(x) ∫ α(x,y) u(y) dy``.
	Thus, this function solves a non-linear system of equations for `u_m` and `u_f`.

	The constraints ``0 ≤ u ≤ ℓ`` are not enforced here, but the outputs of this function
	are truncated in the fixed point iteration loop.
	"""
	function update_singles!(ρ::Float64, δ::Float64, α::Array,
						 ℓ_m::Vector, ℓ_f::Vector, u_m::Vector, u_f::Vector)

		"""
		Vector function representing the system of equations to be solved.
		"""
		function vecsys!(μ::Vector, res::Vector) # vec'd
			# split by sex
			μm = μ[1:length(u_m)]
			μf = μ[length(u_m)+1:end]

			# compute residuals of non-linear system
			mres = δ * ℓ_m - μm .* (δ + ρ * (α * μf))
			fres = δ * ℓ_f - μf .* (δ + ρ * (α' * μm))

			res[:] = [mres; fres] # concatenate into stacked vector
		end # vecsys!

		# initial guess of shares of singles: stacked vector of previous values
		guess = [u_m; u_f]

		# NLsolve
		result = nlsolve(vecsys!, guess)

		u_m[:] = result.zero[1:length(u_m)]
		u_f[:] = result.zero[length(u_m)+1:end]

	end # update_singles!

	"""
		update_values!(ρ, δ, r, α, h, u_m, u_f, w_m, w_f)

	Compute the implied singlehood value functions given a matching function
	and singles distributions.

	For males, the value function equation is:
	``∀x, w(x) = \frac{ρ}{2(r+δ)} ∫ (h(x,y) - w(x) - w(y)) α(x,y) u(y) dy``.
	Thus, this function solves a non-linear system of equations for `w_m` and `w_f`.
	"""
	function update_values!(ρ::Float64, δ::Float64, r::Float64, α::Array, h::Array,
						 u_m::Vector, u_f::Vector, w_m::Array, w_f::Array)
		θ = ρ / (2*(r+δ)) # precompute constant

		"""
		Vector function representing the system of equations to be solved.
		"""
		function vecsys!(ω::Vector, res::Vector) # vec'd
			# split by sex
			ωm = ω[1:length(w_m)]
			ωf = ω[length(w_m)+1:end]

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

		w_m[:] = result.zero[1:length(w_m)]
		w_f[:] = result.zero[length(w_m)+1:end]

	end # update_values!

	"""
	Calculate matching function ``α(x,y)`` from ``S ≥ 0`` condition.

	When `G` is degenerate, this yields the non-random case.
	"""
	function update_match(G::Function, h::Array, w_m::Array, w_f::Array)
		return 1.0 .- G.(-match_surplus(h, w_m, w_f))
	end


	### Helper functions ###

	"""
	Construct production array from function.
	"""
	function prod_array(mtypes::Vector, ftypes::Vector, prodfn::Function)
		h = Array{Float64}(length(mtypes), length(ftypes))

		for (i,x) in enumerate(mtypes), (j,y) in enumerate(ftypes)
			h[i,j] = prodfn(x, y)
		end

		return h

	end # prod_array

	"""
	Construct match surplus array from value functions.
	"""
	function match_surplus(h::Array, w_m::Vector, w_f::Vector)
		S = similar(h)

		for i in 1:length(w_m), j in 1:length(w_f)
			S[i,j] = h[i,j] - w_m[i] - w_f[j]
		end

		return S
	end # match_surplus

end # type
