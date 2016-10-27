using NLsolve
using QuantEcon: compute_fixed_point

type ShimerSmith

	### Parameters ###

	ρ::Float64 # arrival rate of meetings
	δ::Float64 # arrival rate of separations
	r::Float64 # discount rate


	### Exogenous objects ###

	# population distributions
	ℓ_m::Vector
	ℓ_f::Vector

	# production function as array
	h::Array


	### Endogenous equilibrium objects ###

	# masses of singles
	u_m::Vector
	u_f::Vector

	# match function
	α::Array

	# value functions of singles
	V_m::Function
	V_f::Function

	# marital surplus
	S::Array


	### Inner Constructor ###

	"""
		ShimerSmith(ρ, δ, r, ℓ_m, ℓ_f, h)

	Construct a Shimer \& Smith (2000) marriage market model and solve for the equilibrium.
	The equilibrium is the solution to a three-part fixed-point mapping.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64, ℓ_m::Vector, ℓ_f::Vector, h::Array)

		# CHECK: masses of m/f must be valid distro
		if minimum(ℓ_m) < 0.0 || minimum(ℓ_f) < 0.0
			error("Invalid type distribution.")
		end

		# CHECK: masses of m/f match the number of types
		if length(ℓ_m) != size(h)[1] || length(ℓ_f) != size(h)[2]
			error("Number of types inconsistent with production array.")
		end

		# TODO: compute equilibrium objects via fixed point iteration

		# construct instance
		new(ρ, δ, r, ℓ_m, ℓ_f, h, u_m, u_f, α, V_m, V_f, S)

	end # constructor


	### Outer Constructors ###

	"""
		ShimerSmith(ρ, δ, r, Θ_m, Θ_f, ℓ_m, ℓ_f, g)

	Construct a Shimer & Smith (2000) marriage market model using a production function `h(x,y)`.
	"""
	function ShimerSmith(ρ::Float64, δ::Float64, r::Float64,
					  Θ_m::Vector, Θ_f::Vector, ℓ_m::Vector, ℓ_f::Vector, g::Function)
		h = prod_array(Θ_m, Θ_f, g)
		return ShimerSmith(ρ, δ, r, ℓ_m, ℓ_f, h)
	end


	### Equilibrium computation ###

	"""
		update_singles(ρ, δ, ℓ_m, ℓ_f, u_m, u_f, α)

	Compute the implied singles distributions given a matching function from the
	steady-state equilibrium conditions.

	For males, this equation is:
	``∀x, δ(ℓ(x) - u(x)) = ρ u(x) ∫α(x,y)u(y)dy``.
	Thus, this function solves a non-linear system of equations for `u_m` and `u_f`.
	"""
	function update_singles!(ρ::Float64, δ::Float64, α::Array,
						  ℓ_m::Vector, ℓ_f::Vector, u_m::Vector, u_f::Vector)

		"""
		Split vector into male/female, then calculate residuals and recombine.
		"""
		function vecsys!(μ::Vector, res::Vector) # vec'd
			# split by sex
			μm = μ[1:length(u_m)]
			μf = μ[length(u_m)+1:end]

			mres = δ * ℓ_m - μm .* (δ + ρ * (α * μf))
			fres = δ * ℓ_f - μf .* (δ + ρ * (α' * μm))
			res[:] = [vec(mres); vec(fres)] # concatenate into big vector
		end # vecsys!

		# initial guess of shares of singles: stacked vector of previous values
		guess = [u_m; u_f]

		# NLsolve
		result = nlsolve(vecsys!, guess, ftol=1e-16)

		u_m[:] = result.zero[1:length(u_m)]
		u_f[:] = result.zero[length(u_m)+1:end]

		return u_m, u_f
	end # update_singles!

	"""
	Calculate singlehood value functions.
	"""
	function update_values!(ρ::Float64, δ::Float64, r::Float64, α::Array, h::Array,
						 u_m::Vector, u_f::Vector, w_m::Array, w_f::Array)
		# similar system of equations as in update_singles
	end

	"""
	Calculate matching function α from S ≧ 0 condition.
	"""
	function update_match!(h::Array, w_m::Array, w_f::Array, α::Array)
		for i in 1:length(w_m), j in 1:length(w_f)
			α[i,j] = 1*(h[i,j] - w_m[i] - w_f[j] ≥ 0)
		end

		return α
	end

	# TODO: fixed point composition mapping


	### Helper functions ###

	"""
	Construct production array from function.
	"""
	function prod_array(mtypes::Vector, ftypes::Vector, prodfn::Function)
		# generate surp from prodfn and types
		surp = Array{Float64}(length(mtypes), length(ftypes))

		for (i,x) in enumerate(mtypes), (j,y) in enumerate(ftypes)
			surp[i,j] = prodfn(x, y)
		end

		return surp

	end # prod_array

end # type

