include("personalised_k6HFMM_func.jl")

# generate results
simulates, estimates, residuals, glu_bands, ins_bands = gen_results(;test="OGTT");

# save parameter estimates, simulation and residuals
mmts = stack(estimates, [:k1_l, :k1, :k1_r, :k5_l, :k5, :k5_r]);
CSV.write("MMTestimates_2param.csv", mmts; transform=(col, val) -> something(val, missing))
CSV.write("MMTsimulates_2param.csv", simulates)
CSV.write("MMTresiduals_2param.csv", residuals)

# Analyze PLA results
pla_results = CSV.read("pla_results_2param.csv", DataFrame)

# Count unidentifiable parameters
unidentifiable_k1 = count(pla_results.Parameter .== "k1" .&& pla_results.Identifiability .!= "Identifiable")
unidentifiable_k5 = count(pla_results.Parameter .== "k5" .&& pla_results.Identifiability .!= "Identifiable")

println("Unidentifiable k1 parameters: $unidentifiable_k1")
println("Unidentifiable k5 parameters: $unidentifiable_k5")

# Calculate percentage of unidentifiable parameters
total_k1 = count(pla_results.Parameter .== "k1")
total_k5 = count(pla_results.Parameter .== "k5")

percent_unidentifiable_k1 = (unidentifiable_k1 / total_k1) * 100
percent_unidentifiable_k5 = (unidentifiable_k5 / total_k5) * 100

println("Percentage of unidentifiable k1 parameters: $(round(percent_unidentifiable_k1, digits=2))%")
println("Percentage of unidentifiable k5 parameters: $(round(percent_unidentifiable_k5, digits=2))%")