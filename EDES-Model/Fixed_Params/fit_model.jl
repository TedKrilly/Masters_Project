include("func.jl")

try
    # generate results
    OGTTsimulates, OGTTestimates, OGTTresiduals, OGTTglu_bands, OGTTins_bands = gen_results(;test="OGTT")

    # save parameter estimates, simulation and residuals
    ogtts = stack(OGTTestimates, 3:11)
    CSV.write("OGTTestimates.csv", ogtts; transform=(col, val) -> something(val, missing))
    CSV.write("OGTTsimulates.csv", OGTTsimulates)
    CSV.write("OGTTresiduals.csv", OGTTresiduals)

    println("Processing completed successfully!")
catch e
    println("An error occurred:")
    println(e)
    println(stacktrace(catch_backtrace()))
end

