using NLsolve

"""
	StaticMatch(mtypes, ftypes, mdist, fdist, surplus)

Construct a Choo & Siow (2006) marriage market model and solve for the equilibrium.
"""
immutable StaticMatch # for compatibility with Julia 0.5; from 0.6 the keyword is `struct`

	# m/f types
	"vector of vectors of values for each of N male traits"
	mtypes::Vector{Vector}
	"vector of vectors of values for each of L female traits"
	ftypes::Vector{Vector}

	# m/f masses
	"male masses: I_1 x ... x I_N"
	mdist::Array
	"female masses: J_1 x ... x J_L"
	fdist::Array

	"surplus array: (I_1 x ... x I_N) x (J_1 x ... x J_L)"
	surplus::Array

	# equilibrium
	"masses of single males: I_1 x ... x I_N"
	msingle::Array
	"masses of single females: J_1 x ... x J_L"
	fsingle::Array
	"masses of married couples: (I_1 x ... x I_N) x (J_1 x ... x J_L)"
	matches::Array
	"wife's share of marital surplus: (I_1 x ... x I_N) x (J_1 x ... x J_L)"
	wifeshare::Array

	"Inner constructor solves equilibrium and performs sanity checks"
	function StaticMatch(mtypes::Vector{Vector}, ftypes::Vector{Vector},
	                     mdist::Array, fdist::Array, surplus::Array)

		# CHECK: masses of m/f must be proper probability distro
		if minimum(mdist) < 0.0 || minimum(fdist) < 0.0
			error("invalid type distribution")
		end

		# CHECK: masses of m/f match the number of types
		if ndims(mdist) != length(mtypes) || ndims(fdist) != length(ftypes)
			error("type distributions inconsistent with type vectors")
		end

		# compute equilibrium
		msingle, fsingle = equilibrium(surplus, mdist, fdist)
		matches = match_matrix(surplus, msingle, fsingle)
		wifeshare = surplusdiv(matches, msingle, fsingle) ./ (2.0 * surplus)

		# TEST: masses of every match outcome must be strictly positive
		if minimum(msingle) < -1e-5 || minimum(fsingle) < -1e-5
			error("Non-positive mass of singles.")
		end
		if minimum(matches) < -1e-5
			error("Non-positive match mass.")
		end

		# numerical solution may have tiny negatives, set them to zero
		msingle[msingle .< 0.0] = 0.0
		fsingle[fsingle .< 0.0] = 0.0
		matches[matches .< 0.0] = 0.0

		# create instance
		new(mtypes, ftypes, mdist, fdist, surplus, msingle, fsingle, matches, wifeshare)
	end # constructor


end # immutable


"Outer constructor that takes the production function to build the surplus array"
function StaticMatch(men::Vector{Vector}, wom::Vector{Vector},
					 mmass::Array, fmass::Array, prodfn::Function)
	# prodfn(man::Vector, woman::Vector)

	surp = generate_surplus(men, wom, mmass, fmass, prodfn)

	# create instance
	return StaticMatch(men, wom, mmass, fmass, surp)
end

"Outer constructor for one dimensional case"
function StaticMatch(men::Vector, wom::Vector,
					 mmass::Array, fmass::Array, prodfn::Function)
	return StaticMatch(Vector[men], Vector[wom], mmass, fmass, prodfn)
end

"Outer constructor for one dimensional males case"
function StaticMatch(men::Vector, wom::Vector{Vector},
					 mmass::Array, fmass::Array, prodfn::Function)
	return StaticMatch(Vector[men], wom, mmass, fmass, prodfn)
end

"Outer constructor for one dimensional females case"
function StaticMatch(men::Vector{Vector}, wom::Vector,
					 mmass::Array, fmass::Array, prodfn::Function)
	return StaticMatch(men, Vector[wom], mmass, fmass, prodfn)
end



# helper functions

"Compute equilibrium shares of singles."
function equilibrium(surpl::Array, mmass::Array, fmass::Array)

	"""
	The zero of this function gives the equilibrium shares of singles,
	which fully determines the match matrix.
	"""
	function redusys(surp::Array, μm0::Array, μf0::Array)
		mres = similar(μm0) # initialize output array of residuals
		for i in CartesianRange(size(μm0)) #men are i, women are j
			mres[i] = mmass[i] - μm0[i] - ( sfrt(μm0[i]) *
			                                sum([exp(surp[i.I..., j.I...]) * sfrt(μf0[j])
		                                         for j in CartesianRange(size(μf0))]))
		end #for

		fres = similar(μf0) # initialize output array of residuals
		for j in CartesianRange(size(μf0)) #women are j, men are i
			fres[j] = fmass[j] - μf0[j] - ( sfrt(μf0[j]) *
			                                sum([exp(surp[i.I..., j.I...]) * sfrt(μm0[i])
		                                         for i in CartesianRange(size(μm0))]))
		end #for

		return mres, fres
	end # redusys

	function sysvec!(μ0::Vector, res::Vector) # vec'd
		# split and reshape mu
		μm00 = reshape(μ0[1:prod(size(mmass))], size(mmass))
		μf00 = reshape(μ0[prod(size(mmass))+1:end], size(fmass))

		mresid, fresid = redusys(surpl, μm00, μf00)
		res[:] = [vec(mresid); vec(fresid)] # concatenate into big vector
		return res
	end

	# initial guess of shares of singles: stacked vector
	guess = 0.15 * [vec(mmass); vec(fmass)]

	# NLsolve
	result = nlsolve(sysvec!, guess, ftol=1e-16, factor=0.1)

	singmen = reshape(result.zero[1:prod(size(mmass))], size(mmass))
	singwom = reshape(result.zero[prod(size(mmass))+1:end], size(fmass))

	return singmen, singwom
end # equilibrium


function match_matrix(surp::Array, singlemen::Array, singlewom::Array)
	matches = similar(surp)
	for i in CartesianRange(size(matches))
		matches[i] = exp(surp[i]) * sfrt(singlemen[i.I[1:ndims(singlemen)]...] *
								   singlewom[i.I[ndims(singlemen)+1:end]...])
	end # for
	return matches
end # match_matrix

"Saferoot function to prevent DomainError in NLsolve."
function sfrt(x::Real)
	if x < 0.0
		return 0.0
	else
		return sqrt(x)
	end
end  # sfrt

"Construct production array from function."
function generate_surplus(mtypes::Vector{Vector}, ftypes::Vector{Vector},
                          mmass::Array, fmass::Array, prodfn::Function)
	# get dimensions
	Dm = [length(v) for v in mtypes]
	Df = [length(v) for v in ftypes]

	# initialize arrays
	surp = Array{Float64}(Dm..., Df...)
	gent = Vector{Float64}(length(mtypes)) # one man's vector of traits
	lady = Vector{Float64}(length(ftypes))

	for coord in CartesianRange(size(surp))
		for trt in 1:length(mtypes) # loop through traits in coord
			gent[trt] = mtypes[trt][coord[trt]]
		end
		for trt in 1:length(ftypes) # loop through traits in coord
			lady[trt] = ftypes[trt][coord[trt+length(mtypes)]]
		end
		surp[coord] = prodfn(gent, lady)
	end

	return surp

end # generate_surplus

"Wife's consumption out of surplus (aggregate share)."
function surplusdiv(matches::Array, sm::Array, sw::Array)
	share = Array{Float64}(size(matches))

	for i in CartesianRange(size(share))
		share[i] = log(matches[i]) - log(sw[i.I[ndims(sm)+1:end]...])
	end # for

	return share
end # surplus
