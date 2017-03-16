using Distributions

ρ, r, δ, σ = 5.0, 0.05, 0.05, 10.0

@everywhere h(x::Real, y::Real) = x*y # unidimensional production function
@everywhere hsup(x::Vector,y::Vector) = dot(x,y)/2 # supermodular production function
@everywhere hsub(x::Vector,y::Vector) = sqrt(sum(x+y)) # submodular production function

function static_unidim(nmen, nwom, prod)
	
	# types: concave to prevent numerical instability
	men = [log(i+1) for i=1:nmen]
	wom = [log(j+1) for j=1:nwom]

	# masses: unit mass of each sex
	menmass = ones(Float64, nmen)/nmen
	wommass = ones(Float64, nwom)/nwom

	job = @spawn StaticMatch(men, wom, menmass, wommass, prod)
	return job
end

function search_uniform(ntypes, mmass, fmass, prod, σ)
	Θ = Vector(linspace(1.0, 2.0, ntypes)) # types

	# uniform population distributions
	lm = (mmass / ntypes) * ones(Float64, ntypes)
	lf = (fmass / ntypes) * ones(Float64, ntypes)

	job = @spawn SearchClosed(ρ, δ, r, σ, Θ, Θ, lm, lf, prod)
	return job
end

### Static model ###

pam_job = static_unidim(5, 5, hsup) # must be symmetric for the symmetry test
nam_job = static_unidim(3, 5, hsub)

# multidimensional types: symmetric case
n1, n2 = 6, 2
# common type vector
symtypes = Vector[[log(1+i) for i=1:n1], [i for i=1:n2]]
# mass vectors: unit mass of each sex
symdist = ones(Float64, n1, n2)/(n1*n2)

symm_multi_job = @spawn StaticMatch(symtypes, symtypes, symdist, symdist, hsup)

# multidimensional types: asymmetric case
# type vectors
men2 = Vector[[1.0, 1.2, 1.3], [0.0, 1.0]]
wom2 = Vector[[1.0, 1.2, 1.3, 1.35, 1.4], [0.0, 1.0]]
# mass vectors: unit mass of each sex
mdist2 = ones(Float64, 3, 2)/6
fdist2 = ones(Float64, 5, 2)/10

asym_multi_job = @spawn StaticMatch(men2, wom2, mdist2, fdist2, hsup)


### Search model - non-stochastic ###

# instantiate on a worker process
symm_job = search_uniform(20, 100, 100, h, 0)

### Search model - stochastic ###

# symmetric case
rsym_job = search_uniform(20, 100, 100, h, σ)

# asymmetric case: needs σ >~ 10 to converge
rasym_job = search_uniform(20, 50, 100, h, σ)

# multidimensional types: symmetric case
k1, k2 = 10, 2
# common type vector
srsymtypes = Vector[[log(1+i) for i=1:k1], [i for i=1:k2]]
# mass vectors: unit mass of each sex
srsymdist= ones(Float64, k1, k2)/(k1*k2)

srsymm_multi_job = @spawn SearchClosed(ρ, δ, r, σ, srsymtypes, srsymtypes, srsymdist, srsymdist, hsup)
