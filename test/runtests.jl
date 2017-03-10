using MarriageMarkets
using FactCheck

include("setup_tests.jl") # setup common variables

facts("Static model:") do
	context("Unidimensional case") do
		include("test_static_unidim.jl")
	end

	context("Multidimensional case") do
		include("test_static_multidim.jl")
	end
end #facts

facts("Search model:") do
	context("Basic Shimer-Smith model") do
		include("test_search_base.jl")
	end

	context("Endogenous divorce model") do
		include("test_search_stoch.jl")
	end
end #facts

exitstatus() # for CI, makes this code fail if a test fails
