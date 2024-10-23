include("streamlined_func_noins_ident.jl")

# generate results (only for CGM now)
CGMsimulates, CGMestimates, CGMresiduals, CGMglu_bands, CGMins_bands = gen_results(;test="CGM");

# save parameter estimates, simulation and residuals
if !isnothing(CGMsimulates)
    cgms = stack(CGMestimates, 3:14);
    CSV.write("CGMestimates_noins.csv", cgms; transform=(col, val) -> something(val, missing))
    CSV.write("CGMsimulates_noins.csv", CGMsimulates)
    CSV.write("CGMresiduals_noins.csv", CGMresiduals)

    # Read and analyze PLA results
    cgm_pla = CSV.read("pla_results_CGM_noins-Copy.csv", DataFrame)

    # Count unidentifiable parameters
    cgm_unidentifiable = count(cgm_pla.Identifiability .!= "Identifiable")

    println("CGM unidentifiable parameters: $cgm_unidentifiable")

    # Summarize results
    println("\nParameter Identifiability Summary:")
    for param in unique(cgm_pla.Parameter)
        param_results = cgm_pla[cgm_pla.Parameter .== param, :]
        identifiable = count(param_results.Identifiability .== "Identifiable")
        pract_unidentifiable = count(param_results.Identifiability .== "Practically unidentifiable")
        struct_unidentifiable = count(param_results.Identifiability .== "Structurally unidentifiable")
        errors = count(param_results.Identifiability .== "Error in calculation")
        
        println("$param: Identifiable: $identifiable, Practically unidentifiable: $pract_unidentifiable, Structurally unidentifiable: $struct_unidentifiable, Errors: $errors")
    end

    # Print detailed error information
    println("\nDetailed Error Information:")
    error_cases = cgm_pla[cgm_pla.Identifiability .== "Error in calculation", :]
    for row in eachrow(error_cases)
        println("ID: $(row.ID), Parameter: $(row.Parameter)")
        println("Error: $(row.Error_Info)")
        println()
    end
else
    println("No CGM data to process and save.")
end