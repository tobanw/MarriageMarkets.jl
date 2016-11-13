using MarriageMarkets
using FactCheck

facts("Running unit tests.") do
	context("Static model") do
		include("test_static_unidim.jl")
		include("test_static_multidim.jl")
	end
	context("Search model") do
		include("test_search.jl")
	end
end #facts

exitstatus() # for CI, makes this code fail if a test fails
