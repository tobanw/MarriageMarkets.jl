# symmetric case

n1 = 6
n2 = 2

men = Vector[[log(1+i) for i=1:n1], [i for i=1:n2]]
wom = Vector[[log(1+i) for i=1:n1], [i for i=1:n2]]

# mass vectors: unit mass of each sex
mdist = ones(Float64, n1, n2)/(n1*n2)
fdist = ones(Float64, n1, n2)/(n1*n2)

h(x,y) = dot(x,y)/2 #production function

mgmkt = MarriageMarket(men, wom, mdist, fdist, h)

# check symmetry
@fact mgmkt.msingle --> roughly(mgmkt.fsingle)


# asymmetric case

men2 = Vector[[1.0, 1.2, 1.3], [0.0, 1.0]]
wom2 = Vector[[1.0, 1.2, 1.3, 1.35, 1.4], [0.0, 1.0]]

# mass vectors: unit mass of each sex
mdist2 = ones(Float64, 3, 2)/6
fdist2 = ones(Float64, 5, 2)/10

mgmkt2 = MarriageMarket(men2, wom2, mdist2, fdist2, h)

# check that supermodular surplus results in positive assortative matching
@fact mgmkt2.matches[1,1,1,1] + mgmkt2.matches[end,end,end,end] --> greater_than(mgmkt2.matches[1,1,end,end] + mgmkt2.matches[end,end,1,1])
