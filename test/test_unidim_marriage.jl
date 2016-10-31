function basicunidim(nmen, nwom, prod)
	
	# types: concave to prevent numerical instability
	men = [log(i+1) for i=1:nmen]
	wom = [log(j+1) for j=1:nwom]

	# masses: unit mass of each sex
	menmass = ones(Float64, nmen)/nmen
	wommass = ones(Float64, nwom)/nwom

	return ChooSiow(men, wom, menmass, wommass, prod)
end


# instantiate a ChooSiow and verify that super/sub-modularity of production implies positive/negative assortative matching.
n = 5
h(x,y) = dot(x,y)/2 # supermodular production function
pam = basicunidim(n, n, h) # must be symmetric for the symmetry test
nmen = 3
nwom = 5
g(x,y) = sqrt(sum(x+y)) # submodular production function
nam = basicunidim(nmen, nwom, g)

# sum of diag corners greater than sum of anti-diag corners
@fact pam.matches[1,1] + pam.matches[end,end] --> greater_than(pam.matches[1,end] + pam.matches[end,1])
@fact nam.matches[1,1] + nam.matches[end,end] --> less_than(nam.matches[1,end] + nam.matches[end,1])

# symmetry
@fact pam.matches --> roughly(pam.matches')
@fact pam.msingle --> roughly(pam.fsingle)
