using NLsolve
using QuantEcon: compute_fixed_point
using Distributions

const STDNORMAL = Normal()

"""
	SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h; verbose=false, step=0.2)

Construct a Shimer & Smith (2000) marriage market model and solve for the equilibrium.
When match-specific shocks are included, the divorce process is endogenized as in Goussé (2014).

The equilibrium is the solution to a nested fixed-point mapping.
The inner fixed-point solves for a matching equilibrium: consistent strategies given fixed
singles densities.
The outer fixed-point solves for a market equilibrium: singles densities that are consistent
with strategies.

Model selection depends which arguments are provided:
* Match-specific additive shocks ``z ~ N(0, σ)`` when ``σ > 0``
	* Note: some randomness is required for the fixed point iteration to converge with sex asymmetry
* Closed system or inflow/outflow:
	* `ℓ_m, ℓ_f` exogenous: population circulates between singlehood and marriage, no birth/death
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
	v_m::Vector
	"female singlehood (average) value function as vector"
	v_f::Vector

	"match function as array"
	α::Array

	"marital (average) surplus function as array"
	s::Array


	### Inner Constructor ###

	"""
	Inner constructor that accomodates several model variants, depending on the arguments provided.
	It is not meant to be called directly -- instead, outer constructors should call this
	constructor with a full set of arguments, using zero values for unwanted components.
	"""
	function SearchMatch(ρ::Real, δ::Real, r::Real, σ::Real,
	                     γ_m::Vector, γ_f::Vector, ψ_m::Vector, ψ_f::Vector,
	                     ℓ_m::Vector, ℓ_f::Vector, h::Array;
	                     verbose=false, step=0.2)

		### Model Selection ###

		N_m = length(γ_m)
		N_f = length(γ_f)

		if sum(ψ_m) > 0 && sum(ψ_f) > 0 # inflow/outflow model: if death rates provided
			INFLOW = true
		else
			INFLOW = false
		end

		if σ == 0
			STOCH = false
			step = 1 # don't shrink fixed point iteration steps
		else
			STOCH = true
		end


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
		elseif size(h)[1] != N_m || size(h)[2] != N_f
			error("Number of types inconsistent with production array.")
		elseif length(ℓ_m) ≠ N_m || length(ℓ_f) ≠ N_f
			error("Inconsistent number of types.")
		end


		### Compute Constants and Define Functions ###

		if INFLOW
			# population size given directly by inflow and outflow rates
			ℓ_m = γ_m ./ ψ_m
			ℓ_f = γ_f ./ ψ_f
		end

		"cdf of match-specific marital productivity shocks"
		G(x::Real) = STOCH ? cdf(Normal(0, σ), x) : Float64(x ≥ 0) # bool as float

		"μ function, using inverse Mills ratio: E[z|z>q] = σ ϕ(q/σ) / (1 - Φ(q/σ))."
		function μ(a::Real)
			st = quantile(STDNORMAL, 1-a) # pre-compute -s/σ = Φ^{-1}(1-a)
			return σ * (pdf(STDNORMAL, st) - a * st)
		end


		### Steady-State Equilibrium Conditions ###
		"""
		Update population distributions.

		Compute the implied singles distributions given a matching function.
		The steady state conditions equating flows into and out of marriage
		in the deterministic setting are
		```math
		∀x, (δ + ψ(x))(ℓ(x) - u(x)) = ρ u(x) ∫ α(x,y) u(y) dy
		```
		and with productivity shocks and endogenous divorce they are
		```math
		∀x, ℓ(x) - u(x) = ρ u(x) ∫\frac{α(x,y)}{ψ(x)+ψ(y)+δ(1-α(x,y))}u(y)dy
		```
		Thus, this function solves a non-linear system of equations for `u_m` and `u_f`.
		The constraints ``0 ≤ u ≤ ℓ`` are not enforced here, but the outputs of this function
		are truncated in the fixed point iteration loop.
		"""
		function steadystate!(u::Vector, res::Vector) # stacked vector
			# uses the overwritable α in the outer scope
			um, uf = sex_split(u, N_m)

			if STOCH
				αQ = α ./ (δ*(1-α) .+ ψ_m .+ ψ_f') # precompute the fixed integration weights

				# compute residuals of non-linear system
				mres = ℓ_m - um .* (1 .+ ρ .* (αQ * uf))
				fres = ℓ_f - uf .* (1 .+ ρ .* (αQ' * um))
			else # deterministic case
				mres = (δ .+ ψ_m) .* ℓ_m - um .* ((δ .+ ψ_m) .+ ρ .* (α * uf))
				fres = (δ .+ ψ_f) .* ℓ_f - uf .* ((δ .+ ψ_f) .+ ρ .* (α' * um))
			end

			res[:] = [mres; fres] # concatenate into stacked vector
		end # steadystate!

		"""
		Update singlehood value functions for deterministic case only.

		Compute the implied singlehood value functions given a matching function
		and singles distributions.
		```math
		∀x, 2v(x) = ρ ∫ α(x,y) S(x,y) u(y) dy``,\\
		S(x,y) = \frac{h(x,y) - v(x) - v(y)}{r+δ+ψ_m(x)+ψ_f(y)}
		```
		This function solves a non-linear system of equations for the average value
		functions, `v(x) = (r+ψ(x))V(x)`.
		"""
		function valuefunc_base!(v::Vector, res::Vector, u_m::Vector, u_f::Vector, A::Array)
			vm, vf = sex_split(v, N_m)

			# precompute the fixed weights
			αS = A .* match_surplus(vm, vf, A) ./ (r + δ + ψ_m .+ ψ_f')

			# compute residuals of non-linear system
			mres = 2*vm - ρ * (αS * u_f)
			fres = 2*vf - ρ * (αS' * u_m)

			res[:] = [mres; fres] # concatenate into stacked vector
		end # valuefunc_base!

		"""
		Update matching function ``α(x,y)`` from ``S ≥ 0`` condition.

		When `G` is degenerate, this yields the non-random case.
		"""
		function update_match(v_m::Array, v_f::Array, A::Array)
			return 1 .- G.(-match_surplus(v_m, v_f, A)) # in deterministic case, G is indicator function
		end

		"Compute average match surplus array ``s(x,y)`` from value functions."
		function match_surplus(v_m::Vector, v_f::Vector, A::Array)
			if STOCH
				s = h .- v_m .- v_f' .+ δ * μ.(A) ./ (r + δ + ψ_m .+ ψ_f')
			else # deterministic Shimer-Smith model
				s = h .- v_m .- v_f'
			end
			return s
		end


		### Equilibrium Solver ###

		# Initialize guesses for v (deterministic case only): overwritten and reused
		#   in the inner fixed point iteration.
		v_m = 0.5 * h[:,1]
		v_f = 0.5 * h[1,:]

		# Initialize matching array: overwritten and reused in the outer fixed point iteration
		α = 0.5 * ones(Float64, N_m, N_f)

		# rough initial guess for positive assortativity
		#for i in 1:N_m, j in 1:N_f
		#	α[i,j] = 1/(1 + exp(abs(i-j)/10))
		#end

		"""
		Matching equilibrium fixed point operator ``T_α(α)``.
		Solves value functions and match probabilities given singles distributions.
		The value functional equations are:
		```math
		∀x, (r+ψ(x))V(x) = ρ/2 ∫\frac{μ(α(x,y))}{r+δ+ψ(x)+ψ(y)}u(y)dy,\\
		μ(α(x,y)) = α(x,y)(-G^{-1}(1-α(x,y)) + E[z|z > G^{-1}(1-α(x,y))])
		```
		Alternatively, this could be written as an array of point-wise equations
		and solved for α.

		They keyword argument `step` controls the step size of the fixed point iteration.
		Steps must be shrunk or else the iterates can get stuck in an oscillating pattern.
		"""
		function fp_matching_eqm(A::Array, u_m::Vector, u_f::Vector)
			# overwrite v_m, v_f to reuse as initial guess for nlsolve
			if STOCH
				μα = μ.(A) ./ (r + δ + ψ_m .+ ψ_f') # precompute μ/deno

				v_m[:] = 0.5*ρ * (μα * u_f)
				v_f[:] = 0.5*ρ * (μα' * u_m)

			else # need to solve non-linear system because α no longer encodes s
				v_m[:], v_f[:] = sex_solve((x,res)->valuefunc_base!(x, res, u_m, u_f, A), v_m, v_f)
			end

			# shrink update step size
			return (1-step) .* A .+ step * update_match(v_m, v_f, A)
		end

		"""
		Market equilibrium fixed point operator ``T_u(u_m, u_f)``.
		Solves for steady state equilibrium singles distributions, with strategies
		forming a matching equilibrium.
		"""
		function fp_market_eqm(u::Vector; inner_tol=1e-4)

			um, uf = sex_split(u, N_m)

			# nested fixed point: overwrite α to be reused as initial guess for next call
			α[:] = compute_fixed_point(x->fp_matching_eqm(x, um, uf), α,
			                           err_tol=inner_tol, verbose=verbose)

			# steady state distributions
			# TODO: precompute αQ in outer scope?
			um_new, uf_new = sex_solve(steadystate!, um, uf)

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

		# Solve fixed point

		# fast rough compututation of equilibrium by fixed point iteration
		u_fp0 = compute_fixed_point(fp_market_eqm, 0.1*[ℓ_m; ℓ_f],
		                            print_skip=10, verbose=1+verbose) # initial guess u = 0.1*ℓ

		um_0, uf_0 = sex_split(u_fp0, N_m)

		# touch up with high precision fixed point solution
		u_fp = compute_fixed_point(x->fp_market_eqm(x, inner_tol=1e-8), [um_0; uf_0],
		                           err_tol=1e-8, verbose=1+verbose)
		
		u_m, u_f = sex_split(u_fp, N_m)

		s = match_surplus(v_m, v_f, α)

		# construct instance
		new(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, h, ℓ_m, ℓ_f, u_m, u_f, v_m, v_f, α, s)

	end # constructor

end # type


### Outer Constructors ###

# Recall inner constructor: SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h)

"""
	SearchClosed(ρ, δ, r, σ, Θ_m, Θ_f, ℓ_m, ℓ_f, g; verbose=false, step=0.2)

Constructs marriage market equilibrium of closed-system model with match-specific productivity shocks and production function ``g(x,y)``.
"""
function SearchClosed(ρ::Real, δ::Real, r::Real, σ::Real,
                      Θ_m::Vector, Θ_f::Vector, ℓ_m::Vector, ℓ_f::Vector,
                      g::Function; verbose=false, step=0.2)
	# irrelevant arguments to pass as zeros
	ψ_m = zeros(ℓ_m)
	ψ_f = zeros(ℓ_f)
	γ_m = zeros(ℓ_m)
	γ_f = zeros(ℓ_f)

	h = prod_array(Θ_m, Θ_f, g)
	return SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h; verbose=verbose, step=step)
end

"""
	SearchInflow(ρ, δ, r, σ, Θ_m, Θ_f, γ_m, γ_f, ψ_m, ψ_f, g; verbose=false, step=0.2)

Constructs marriage market equilibrium of inflow and death model with match-specific productivity shocks and production function ``g(x,y)``.
"""
function SearchInflow(ρ::Real, δ::Real, r::Real, σ::Real,
                      Θ_m::Vector, Θ_f::Vector, γ_m::Vector, γ_f::Vector,
                      ψ_m::Vector, ψ_f::Vector, g::Function;
					  verbose=false, step=0.2)
	# irrelevant arguments to pass as zeros
	ℓ_m = zeros(γ_m)
	ℓ_f = zeros(γ_f)

	h = prod_array(Θ_m, Θ_f, g)
	return SearchMatch(ρ, δ, r, σ, γ_m, γ_f, ψ_m, ψ_f, ℓ_m, ℓ_f, h; verbose=verbose, step=step)
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
