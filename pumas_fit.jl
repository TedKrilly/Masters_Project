using Pumas
using CSV
using DataFrames

# Load model
include("pumas_func.jl")

# Load data
data = CSV.read("data_file.csv", DataFrame)

# Create Population object
pop = read_pumas(
    data,
    id = :ID,
    time = :time,
    observations = [:glucose],
    covariates = [:weight, :data_type],
    event_data = false
)

# Set initial parameter values
init_params = (
    k1 = (1.35e-2, 0.5),
    k5 = (3.80e-3, 0.5),
    k6 = (5.82e-1, 0.5),
    τ_g = (2.5, 0.5),
    Ω = Diagonal([0.1, 0.1, 0.1, 0.1]),
    σ = (0.1,)
)

# Fit the model
fit_result = fit(
    edes_model,
    pop,
    init_params,
    SAEM(
        iterations = (1000, 500, 500),
        nchain = 1,
        tol = 1e-3,
    )
)

# Print results
println("Parameter Estimates:")
println(coef(fit_result))

println("\nConfidence Intervals:")
infer_result = infer(fit_result)
println(coeftable(infer_result))