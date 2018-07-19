# Import dependency. 
using BenchmarkTools, Expectations, Distributions

#Create benchmark group and benchmarks
benchmarks = BenchmarkGroup()

#Put in specific benchmarks
dist1 = DiscreteUniform(1, 10)
E1 = expectation(dist1)
benchmarks["discreteuniform"] = @benchmarkable E1($(x -> x))
dist2 = Normal()
E2 = expectation(dist2)
benchmarks["standardnormal"] = @benchmarkable E2($(x -> x))
dist3 = Exponential(2)
E3 = expectation(dist3)
benchmarks["exponential"] = @benchmarkable E3($(x -> x))
#...

results = run(benchmarks) # Get results. 
results = median(results) # Condense to median. 

# To save results, manually call in the REPL: BenchmarkTools.save("benchmarks.json", results)

#Compare to old results
try
  oldresults= BenchmarkTools.load("benchmarks.json")[1]
  judge(oldresults, results)
catch err
  error("Couldn't load file- make sure that you've previously saved results.", err.prefix)
end