using MarriageMarkets
using FactCheck

facts("Running unit tests.") do
	context("Marriage market equilibrium.") do
		include("test_unidim_marriage.jl")
		include("test_multidim_marriage.jl")
	end
end #facts

exitstatus() # for CI, makes this code fail if a test fails
