using FactCheck

# compute models in parallel: multiprocess
addprocs(Sys.CPU_CORES) # add a worker process per core
print_with_color(:white, "Setup:\n")
print("  > Using $(nprocs()-1) worker processes\n")

using MarriageMarkets

print("  > Computing models...\n")
include("setup_tests.jl") # multiprocess model computations

facts("Static model:") do
	context("One-dimensional types") do
		include("test_static_unidim.jl")
	end

	context("Multidimensional types") do
		include("test_static_multidim.jl")
	end
end #facts

facts("Search model:") do
	context("Basic Shimer-Smith model") do
		include("test_search_base.jl")
	end

	context("Stochastic model") do
		include("test_search_stoch.jl")
	end
end #facts

exitstatus() # for CI, makes this code fail if a test fails
