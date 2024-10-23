export edes_delayfix!, loss, set_conditions, update_conditions, PLA, plot_bands, gen_results
using Pkg
using CSV
using DataFrames
using Parameters: @unpack
using ComponentArrays
using DifferentialEquations
using SciMLSensitivity
using UnPack
using Optimization
using LineSearches
using OptimizationNLopt
using OptimizationMultistartOptimization
using ForwardDiff
using LikelihoodProfiler
using Plots

function load_k6_estimates()
    k6_estimates = CSV.read("k6_estimates.csv", DataFrame)
    return Dict(zip(k6_estimates.ID, k6_estimates.value))
end

function edes_delayfix!(du, u, p, t, c_)
    # model equations
    @unpack p_fixed, c, input = c_

    Mg, Gpl, Ipl, Gi = u

    mgmeal = p_fixed[10] * p[1]^p_fixed[10] * t^(p_fixed[10] - 1) * exp(-(p[1] * t)^p_fixed[10]) * input[1]
    mgpl = p_fixed[2] * Mg

    # glucose in gut
    du[1] = mgmeal - mgpl

    c11 = 0.043 * (p_fixed[11] + p_fixed[12]) / p_fixed[12]
    du[4] = (1.0/2.5) * (Gpl - Gi)
    gliv = c[3] - p_fixed[3] * (Gpl - p_fixed[12]) - p_fixed[4] * c[4] * (Ipl - p_fixed[13])
    ggut = p_fixed[2] * (c[1] / (c[2] * input[2])) * Mg
    gnonit = c11 * (Gpl / (p_fixed[11] + Gpl)) 
    git = p[2] * c[4] * Ipl * (Gpl / (p_fixed[11] + Gpl))
    gren = abs(((c[10] / (c[2] * input[2]) * (Gpl - c[8]))) * tanh(100 * (Gpl - c[8])) / 2 + (c[10] / (c[2] * input[2]) * (Gpl - c[8])) / 2)
   
    # glucose in plasma    
    du[2] = gliv + ggut - gnonit - git - gren

    ipnc = (c[4]^-1) * (p_fixed[6] * (Gpl - p_fixed[12]) + (p_fixed[7] / c[5]) * p_fixed[12] + (p_fixed[8] * c[6]) * du[2])
    iliv = c[12] * Ipl
    iif = p_fixed[9] * (Ipl - p_fixed[13])

    # insulin in plasma
    du[3] = ipnc - iliv - iif
    nothing
end;

function loss(p, prob, data)
    # loss function

    # load measurements and their timestamps
    data_glu = data[1,findall(!ismissing, data[1,4:end])]

    ts_glu = parse.(Float64, names(data_glu))

    prob = remake(prob)
    sol = solve(prob, RadauIIA5(), p=p, saveat=1; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
    
    loss = 0.0
    if any((s.retcode != :Success for s in sol))
        loss = Inf
    else
        loss_g = sum(abs2, ((stack(sol(ts_glu).u)[2,:] - Vector(data_glu))/maximum(Vector(data_glu))))
    
        loss = 100*loss_g
    end
    return loss
end;

function set_conditions(test, k6)
    # parameters to estimate & upper and lower bounds for parameter estimation
    p_adj = [
        1.35e-2,  #k1 #0.01
        3.80e-3,  #k5 #0.04
        #k6 removed from here
    ]

    lb = [1e-7, 1e-7]
    ub = [0.3, 1.0]
    
    # fixed parameters/constants
    p_fixed = [
        1.35e-2,  #k1
        6.33e-1,  #k2
        5.00e-5,  #k3
        1.00e-3,  #k4
        3.80e-3,  #k5
        k6,  #k6 is now a variable  that is created from the laod_k6_estimates function and used in the gen_results function
        1.15,     #k7
        4.71,     #k8
        1.08e-2,  #k9
        1.35,     #sigma /index 10
        0.63,     #Km /index 11
        5.0,      #Gb /index 12
        8.0,      #Ib /index 13
        5.0,      #gpllowerbound
        2.84e-1
    ]

    c = [
        0.005551, #f index 1
        17 / 70, # vg index 2
        0.043, # gbliv index 3
        1.0, # beta index 4
        31.0, # taui index 5
        3.0, # taud index 6
        13.0 / 70.0, # vi index 7
        9.0, # Gthpl index 8
        30.0, # t_integralwindow index 9
        0.1, # c1 index 10
        0.043 * (p_fixed[11] + p_fixed[12]) / p_fixed[12] - p_adj[2] * 1 * p_fixed[13],
        p_fixed[7] * (p_fixed[12] / (31 * p_fixed[13])) 
    ]

    input = [
        75000.0, # D
        70.0, # Mb
        0.0 # t_meal_start
    ]

    p_fix = ComponentArray(p_fixed=p_fixed, c=c, input=input)
    
    # simulation time span
    tspan = (0.0, 240.0)    
    u0 = [0.0, p_fixed[12], p_fixed[13], p_fixed[12]]

    return p_adj, lb, ub, p_fix, u0, tspan
end;

function update_conditions!(p_adj, p_fix, u0, data; test="OGTT")
    # used to update the conditions in set_condition to the particular data used in calibration

    u0[2] = data[1, 4] # Glucose, tp=0
    u0[3] = data[2, 4] # Insulin, tp=0
    u0[4] = data[1, 4] # interstitial glucose tp=0

    p_fix.p_fixed[12] = data[1, 4]  # update to measured basal glucose
    p_fix.p_fixed[13] = data[1, 4]  # update to measured basal insulin

    return u0
end;

function PLA(opt, lb, ub, prob, data, storage; id, output_file="pla_results_2param.csv")
    scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100)]
    theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200)]
    param_names = ["k1", "k5"]
   
    α = loss(opt, prob, data) + 3.84

    num_params = length(opt)
    intervals_ = Vector{ParamInterval}(undef, num_params)
    results = DataFrame(
        ID = String[],
        Parameter = String[],
        Lower_Bound = Float64[],
        Upper_Bound = Float64[],
        Identifiability = String[]
    )

    for i in 1:num_params
        param_name = param_names[i]
        try
            intervals_[i] = get_interval(
                opt,
                i,
                x->loss(x,prob,data),
                :CICO_ONE_PASS,
                loss_crit=α,
                theta_bounds=theta_bounds,
                scan_bounds=scan_bounds[i],
                scale=fill(:direct, length(opt)),
                scan_tol=1e-5,
                local_alg = :LN_NELDERMEAD,
                max_iter=10^5
            )
           
            lower_bound = intervals_[i].result[1].value
            upper_bound = intervals_[i].result[2].value
            storage[1, i*3-2] = lower_bound
            storage[1, i*3] = upper_bound
            storage[1, i*3-1] = opt[i]

            param_range = ub[i] - lb[i]
            interval_width = upper_bound - lower_bound
            if interval_width < 0.1 * param_range
                identifiability = "Identifiable"
            elseif interval_width < param_range
                identifiability = "Practically unidentifiable"
            else
                identifiability = "Structurally unidentifiable"
            end
            push!(results, (id, param_name, lower_bound, upper_bound, identifiability))
        catch e
            println("Error in PLA for parameter $param_name of ID $id: $e")
            push!(results, (id, param_name, NaN, NaN, "Error in calculation"))
            storage[1, i*3-2] = NaN
            storage[1, i*3-1] = NaN
            storage[1, i*3] = NaN
        end
    end

    CSV.write(output_file, results, append=true)
    println("PLA results for ID $id saved to $output_file")

    return results
end

function bands(opt, sol, lb, ub, prob, data, ts_sim; test="OGTT")
    # estimate confidence intervals around the observables as a function of parameters over time.

    # simulate observable over time
    function bands_g(p, i) # plasma glu
        prob = remake(prob) 
        sol_ = solve(prob, RadauIIA5(), p=p, saveat=1; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        return sol_[2,i]
    end;

    function bands_i(p, i) # plasma ins
        prob = remake(prob) 
        sol_ = solve(prob, RadauIIA5(), p=p, saveat=1; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        return sol_[3,i]
    end;
    
    # scan and parameter bounds
    scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100)]
    theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200)]

    # confidence threshold
    α = loss(opt, prob, data) + 3.84

    # estimate confidence intervals with LikelihoodProfiler.jl
    intervals_g = Vector{ParamInterval}(undef,length(ts_sim))
    for i in eachindex(ts_sim)
        @time intervals_g[i] = get_interval(
            opt,
            p->bands_g(p,i),
            x->loss(x,prob,data),
            :CICO_ONE_PASS,
            loss_crit = α,
            theta_bounds=theta_bounds,
            scan_bounds=(0.00001,500.0),
            scan_tol = 1e-5,
            local_alg = :LN_NELDERMEAD,
        )
    end

    glbands = [iv.result[1].value for iv in intervals_g]
    gubands = [iv.result[2].value for iv in intervals_g]

    u_g = Array(sol)[2,1:131];

    intervals_i = Vector{ParamInterval}(undef,length(ts_sim))
    for i in eachindex(ts_sim)
        @time intervals_i[i] = get_interval(
            opt,
            p->bands_i(p,i),
            x->loss(x,prob,data),
            :CICO_ONE_PASS,
            loss_crit = α,
            theta_bounds=theta_bounds,
            scan_bounds=(0.00001,500.0),
            scan_tol = 1e-5,
            local_alg = :LN_NELDERMEAD,
        )
    end

    # save confidence bands (confidence interval of observables over time)
    ilbands = [iv.result[1].value for iv in intervals_i]
    iubands = [iv.result[2].value for iv in intervals_i]
    return glbands, gubands, ilbands, iubands
end;

function plot_bands(sim, data, glbands, gubands, ilbands, iubands; test)
    # plot model simulation with confidence intervals

    u_g = Array(sim)[2,1:241];
    u_i = Array(sim)[3,1:241];

    data_glu = data[1,findall(!ismissing, data[1,4:end])]
    data_ins = data[2,findall(!ismissing, data[2,4:end])]   
    ts_glu = parse.(Float64, names(data_glu)) 
    ts_ins = parse.(Float64, names(data_ins))
    
    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 2)
    p[1] = Plots.plot(Float64.(0:240), u_g, xlabel="time", ylabel="glucose, mmol/L", labels="Sim", line=3,
    ribbon = (u_g-glbands, gubands-u_g), fc=:orange, fa=0.3);
    Plots.scatter!(ts_glu, Vector(data_glu), labels="data", legend=:outertopright);
    p[2] = Plots.plot(Float64.(0:240), u_i, xlabel="time", ylabel="insulin, mU/L", labels="Sim", line=3,
    ribbon = (u_i-ilbands, iubands-u_i), fc=:orange, fa=0.3);
    Plots.scatter!(ts_ins, Vector(data_ins), labels="data", legend=:outertopright);
    p = Plots.plot(Plots.plot(p[1],p[2],layout=(2, 1)), legendfontsize=8, legend=:outertopleft);
    display(p)
    return p
end;

function gen_results(;test="OGTT")
    println("Starting gen_results function")
    data = CSV.read("clean_mmt_glucose_and_insulin_data.csv", DataFrame)
    println("Data loaded. Shape: ", size(data))
    
    # select OGTT data
    data = data[data.test .== "OGTT", :]
    sort!(data, :time)
    # convert to wide format
    data_wide = unstack(data, [:ID, :metab, :Condition], :time, :VAL)
    println("Data processed. Shape of data_wide: ", size(data_wide))

    # group for iterating parameter estimation over individuals
    grp_data = groupby(data_wide, [:ID, :Condition])
    println("Number of groups: ", length(grp_data))

    # create DF to store simulated glucose and insulin
    simulated = copy(data_wide[:,1:3])
    sim_ts = DataFrame(missings(Float64,nrow(simulated),241), :auto)
    simulated = hcat(simulated, sim_ts)
    glu_bands = copy(simulated)
    ins_bands = copy(simulated)
    grp_simulated = groupby(simulated, [:ID, :Condition])
    grp_glu_bands = groupby(glu_bands, [:ID, :Condition])
    grp_ins_bands = groupby(ins_bands, [:ID, :Condition])

    # create DF to store parameter estimates
    estimates = copy(unique(data_wide[:,[:ID, :Condition]]))
    n_rows = nrow(estimates)
    println("Number of rows in estimates: ", n_rows)
    est_p = DataFrame(
        k1_l = Vector{Any}(missing, n_rows),
        k1 = Vector{Any}(missing, n_rows),
        k1_r = Vector{Any}(missing, n_rows),
        k5_l = Vector{Any}(missing, n_rows),
        k5 = Vector{Any}(missing, n_rows),
        k5_r = Vector{Any}(missing, n_rows)
    )
    estimates = hcat(estimates, est_p)
    println("Shape of estimates after adding parameters: ", size(estimates))
    grp_estimates = groupby(estimates, [:ID, :Condition])

    # Load k6 estimates
    k6_estimates = load_k6_estimates()
    println("k6 estimates loaded. Number of estimates: ", length(k6_estimates))

    # create DF to store SSR
    residuals = copy(data_wide)
    residuals[:,(4:end)] .= missing
    grp_residuals = groupby(residuals, [:ID, :Condition])
    
    n = length(grp_data)

    for i in 1:n
        println("Processing group $i of $n")
        println(grp_data[i][1,:ID], grp_data[i][1,:Condition])

        # Get the personalized k6 value for this individual
        id = grp_data[i][1,:ID]
        k6 = get(k6_estimates, id, 5.82e-1)  # Use default value if not found
        println("k6 value for ID $id: $k6")
        
        # Set conditions with personalized k6
        p_adj, lb, ub, p_fix, u0, tspan = set_conditions(test, k6)

        u0 = update_conditions!(p_adj, p_fix, u0, grp_data[i]; test)
        prob = ODEProblem((du, u, p, t) -> edes_delayfix!(du, u, p, t, p_fix), u0, tspan, p_adj)
        optf = OptimizationFunction((x, p_adj) -> loss(x, prob, grp_data[i]))
        opt_prob = OptimizationProblem(optf, p_adj, lb=lb, ub=ub)
        opt_sol = solve(opt_prob, MultistartOptimization.TikTak(1000), NLopt.LN_NELDERMEAD(), abstol = 1e-7, reltol = 1e-7)
        println("Optimization done. Solution: ", opt_sol)

        _prob = remake(prob, p=opt_sol.u)
        sol_opt = solve(_prob, RadauIIA5(), saveat=1; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        sol = stack(sol_opt.u)

        # Create storage variable for PLA
        storage = fill(NaN, 1, 6)  # 2 parameters * 3 values each (lower, estimate, upper)

        # Perform PLA
        results = PLA(opt_sol.u, lb, ub, _prob, grp_data[i], storage; 
                      id=id, output_file="pla_results_2param.csv")
        
        # Store PLA results
        grp_estimates[i][1, [:k1_l, :k1, :k1_r, :k5_l, :k5, :k5_r]] = storage

        # store simulated glucose and insulin
        grp_simulated[i][1,4:end] = Array(sol)[2,1:241]
        grp_estimates[i][1, :k1] = opt_sol.u[1]
        grp_estimates[i][1, :k5] = opt_sol.u[2]
        grp_simulated[i][2,4:end] = Array(sol)[3,1:241]

        ts_glu = Int64.(parse.(Float64, names(grp_residuals[i][isequal.(grp_residuals[i].metab,"gluc"), 4:end]))) .+1
        ts_ins = Int64.(parse.(Float64, names(grp_residuals[i][isequal.(grp_residuals[i].metab,"ins"), 4:end]))) .+1

        # store residuals
        grp_residuals[i][1, 4:end] = Array(grp_data[i][isequal.(grp_data[i].metab,"gluc"), 4:end])[1,:] - Array(sol)[2,ts_glu]
        grp_residuals[i][2, 4:end] = Array(grp_data[i][isequal.(grp_data[i].metab,"ins"), 4:end])[1,:] - Array(sol)[3,ts_ins]

        # You can keep or remove the bands calculation depending on whether you need it
        # glbands, gubands, ilbands, iubands = bands(opt_sol.u, sol_opt, lb, ub, _prob, grp_data[i], Float64.(0:240); test)
        # grp_glu_bands[i][1,4:end] = glbands
        # grp_glu_bands[i][2,4:end] = gubands
        # grp_ins_bands[i][1,4:end] = ilbands
        # grp_ins_bands[i][2,4:end] = iubands
        
        # plot_bands(sol, grp_data[i], glbands, gubands, ilbands, iubands; test)
    end

    println("Processing complete. Cleaning up data...")
    dropmissing!(simulated, :x1)
    dropmissing!(glu_bands, :x1)
    dropmissing!(ins_bands, :x1)
    
    println("Returning results")
    return simulated, estimates, residuals, glu_bands, ins_bands
end