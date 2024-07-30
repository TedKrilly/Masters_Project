using CSV
using DataFrames
include("func_copy.jl")
include("test_robustness2.jl")

# Load  OGTT estimates
OGTTestimates = CSV.read("OGTTestimates.csv", DataFrame)
OGTTestimates_wide = unstack(OGTTestimates, [:ID, :Condition], :variable, :value)
OGTTestimates_final = select(OGTTestimates_wide, [:ID, :Condition, :k1, :k5, :k6])

# Load complex meal data
complex_meal_data = CSV.read("mmt_test.csv", DataFrame)

# Filter OGTTestimates_final to only include IDs present in complex_meal_data
matching_ids = intersect(Set(OGTTestimates_final.ID), Set(complex_meal_data.ID))
OGTTestimates_final = filter(row -> row.ID in matching_ids, OGTTestimates_final)

# Test the robustness of the calibrated model on complex meal data
simulated_complex, residuals_complex, metrics_results = test_robustness(OGTTestimates_final, complex_meal_data)

# Save simulated results and residuals for complex meal data
CSV.write("ComplexMealSimulates.csv", simulated_complex)
CSV.write("ComplexMealResiduals.csv", residuals_complex)

# Print metrics for each subject and condition
for (key, metrics) in metrics_results
    println("Metrics for $key:")
    for (metric, value) in metrics
        println("  $metric: $value")
    end
    println()
end

# save these metrics
metrics_df = DataFrame(ID_Condition = String[], Metric = String[], Value = Float64[])
for (key, metrics) in metrics_results
    for (metric, value) in metrics
        push!(metrics_df, (key, metric, value))
    end
end
CSV.write("ComplexMealMetrics.csv", metrics_df)