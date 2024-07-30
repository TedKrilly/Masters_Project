export test_robustness, calculate_metrics

using Statistics, Trapz

function calculate_metrics(observed, simulated)
    # Time in Range
    function calculate_tir(glucose_values, lower_bound=3.9, upper_bound=10.0)
        in_range = (glucose_values .>= lower_bound) .& (glucose_values .<= upper_bound)
        return mean(in_range) * 100
    end

    # Find peak and time to peak
    function find_peak(time, values)
        peak_value, peak_index = findmax(values)
        peak_time = time[peak_index]
        return peak_value, peak_time
    end

    # Area Under the Curve
    function calculate_auc(time, values)
        return trapz(time, values)
    end

    # R-squared
    function r_squared(observed, predicted)
        ss_res = sum((observed .- predicted).^2)
        ss_tot = sum((observed .- mean(observed)).^2)
        return 1 - (ss_res / ss_tot)
    end

    # Calculate metrics
    observed_tir = calculate_tir(observed[observed.metab .== "gluc", :VAL])
    simulated_tir = calculate_tir(simulated[simulated.metab .== "gluc", :VAL])

    obs_glucose = observed[observed.metab .== "gluc", :]
    sim_glucose = simulated[simulated.metab .== "gluc", :]
    obs_insulin = observed[observed.metab .== "ins", :]
    sim_insulin = simulated[simulated.metab .== "ins", :]

    obs_peak_glucose, obs_peak_glucose_time = find_peak(obs_glucose.time, obs_glucose.VAL)
    sim_peak_glucose, sim_peak_glucose_time = find_peak(sim_glucose.time, sim_glucose.VAL)
    obs_peak_insulin, obs_peak_insulin_time = find_peak(obs_insulin.time, obs_insulin.VAL)
    sim_peak_insulin, sim_peak_insulin_time = find_peak(sim_insulin.time, sim_insulin.VAL)

    obs_glucose_auc = calculate_auc(obs_glucose.time, obs_glucose.VAL)
    sim_glucose_auc = calculate_auc(sim_glucose.time, sim_glucose.VAL)
    obs_insulin_auc = calculate_auc(obs_insulin.time, obs_insulin.VAL)
    sim_insulin_auc = calculate_auc(sim_insulin.time, sim_insulin.VAL)

    glucose_r2 = r_squared(obs_glucose.VAL, sim_glucose.VAL)
    insulin_r2 = r_squared(obs_insulin.VAL, sim_insulin.VAL)

    return Dict(
        "observed_tir" => observed_tir,
        "simulated_tir" => simulated_tir,
        "obs_peak_glucose" => obs_peak_glucose,
        "obs_peak_glucose_time" => obs_peak_glucose_time,
        "sim_peak_glucose" => sim_peak_glucose,
        "sim_peak_glucose_time" => sim_peak_glucose_time,
        "obs_peak_insulin" => obs_peak_insulin,
        "obs_peak_insulin_time" => obs_peak_insulin_time,
        "sim_peak_insulin" => sim_peak_insulin,
        "sim_peak_insulin_time" => sim_peak_insulin_time,
        "obs_glucose_auc" => obs_glucose_auc,
        "sim_glucose_auc" => sim_glucose_auc,
        "obs_insulin_auc" => obs_insulin_auc,
        "sim_insulin_auc" => sim_insulin_auc,
        "glucose_r2" => glucose_r2,
        "insulin_r2" => insulin_r2
    )
end

function test_robustness(calibrated_params, complex_meal_data, use_cgm=false)
    simulated_complex = DataFrame(ID = String[], Condition = String[], test = String[], metab = String[], time = Float64[], VAL = Float64[])
    residuals_complex = copy(complex_meal_data)
    metrics_results = Dict()
   
    grp_data = groupby(complex_meal_data, [:ID, :Condition])
   
    for group in grp_data
        id = group[1, :ID]
        condition = group[1, :Condition]
        
        matching_params = calibrated_params[calibrated_params.ID .== id, :]
        if isempty(matching_params)
            error("No matching parameters found for ID: $id")
        end
        
        if use_cgm
            p_adj = Vector(matching_params[1, [:k1, :k5, :k6, :tau_g]])
        else
            p_adj = Vector(matching_params[1, [:k1, :k5, :k6]])
        end
       
        p_fix = set_conditions("OGTT")[4]
        u0 = [0.0, group[group.metab .== "gluc", :VAL][1], group[group.metab .== "ins", :VAL][1], group[group.metab .== "gluc", :VAL][1]]
        tspan = (0.0, maximum(group.time))
       
        prob = if use_cgm
            ODEProblem((du, u, p, t) -> edes_delay!(du, u, p, t, p_fix), u0, tspan, p_adj)
        else
            ODEProblem((du, u, p, t) -> edes_delayfix!(du, u, p, t, p_fix), u0, tspan, p_adj)
        end
       
        sol = solve(prob, RadauIIA5(), saveat=unique(group.time); abstol=1e-7, reltol=1e-5, maxiters=1e7)
       
        for row in eachrow(group)
            t = row.time
            metab = row.metab
            observed = row.VAL
           
            if metab == "gluc"
                simulated = use_cgm ? sol(t)[4] : sol(t)[2]
            elseif metab == "ins"
                simulated = sol(t)[3]
            else
                error("Unknown metabolite type: $metab")
            end
           
            push!(simulated_complex, (id, condition, "HFMM", metab, t, simulated))
           
            residual = observed - simulated
            residuals_complex[residuals_complex.ID .== id .&& residuals_complex.Condition .== condition .&& residuals_complex.metab .== metab .&& residuals_complex.time .== t, :VAL] .= residual
        end

        # Calculate metrics for this group
        group_simulated = simulated_complex[simulated_complex.ID .== id .&& simulated_complex.Condition .== condition, :]
        metrics = calculate_metrics(group, group_simulated)
        metrics_results["$(id)_$(condition)"] = metrics
    end
   
    return simulated_complex, residuals_complex, metrics_results
end
