include("k6HFMM_func.jl")

# generate results
OGTTsimulates, OGTTestimates, OGTTresiduals, OGTTglu_bands, OGTTins_bands = gen_results(;test="OGTT");
# CGMsimulates, CGMestimates, CGMresiduals, CGMglu_bands, CGMins_bands = gen_results(;test="CGM");

# save parameter estimates, simulation and residuals
ogtts = stack(OGTTestimates, 3:4);
CSV.write("HFMM_k6_estimates.csv", ogtts; transform=(col, val) -> something(val, missing))
CSV.write("HFMM_k6_simulates.csv", OGTTsimulates)
CSV.write("HFMM_k6_residuals.csv", OGTTresiduals)

# cgms = stack(CGMestimates, 3:14);
# CSV.write("CGMestimates.csv", cgms; transform=(col, val) -> something(val, missing))
# CSV.write("CGMsimulates.csv", CGMsimulates)
# CSV.write("CGMresiduals.csv", CGMresiduals)
