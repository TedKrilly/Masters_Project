using CSV
using DataFrames
include("func_copy.jl")
include("test_robustness2.jl")
include("data_preprocessing.jl")  # Include the new preprocessing file

# # Load OGTT estimates
# OGTTestimates = CSV.read("OGTTestimates.csv", DataFrame)
# OGTTestimates_wide = unstack(OGTTestimates, [:ID, :Condition], :variable, :value)
# OGTTestimates_final = select(OGTTestimates_wide, [:ID, :Condition, :k1, :k5, :k6])

# Load CGM estimates
CGMestimates = CSV.read("CGMestimates.csv", DataFrame)
CGMestimates_wide = unstack(CGMestimates, [:ID, :Condition], :variable, :value)
CGMestimates_final = select(CGMestimates_wide, [:ID, :Condition, :k1, :k5, :k6, :tauG])

# Load complex meal data
complex_meal_data = CSV.read("mmt_glucose_insulin_data.csv", DataFrame)

# Preprocess the data
processed_params, processed_meal_data = preprocess_data(OGTTestimates_final, complex_meal_data)

# Print diagnostic information
println("Original OGTT estimates shape: $(size(OGTTestimates_final))")
println("Original complex meal data shape: $(size(complex_meal_data))")
println("Processed OGTT estimates shape: $(size(processed_params))")
println("Processed complex meal data shape: $(size(processed_meal_data))")

# Test the robustness of the calibrated model on processed complex meal data
simulated_complex, residuals_complex, metrics_results = test_robustness(processed_params, processed_meal_data)

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
