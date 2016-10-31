using MarriageMarkets
using FactCheck

facts("Running unit tests.") do
	context("Choo & Siow model") do
		include("test_unidim_marriage.jl")
		include("test_multidim_marriage.jl")
	end
	context("Shimer & Smith model") do
		include("test_shimer-smith.jl")
	end
end #facts

exitstatus() # for CI, makes this code fail if a test fails
