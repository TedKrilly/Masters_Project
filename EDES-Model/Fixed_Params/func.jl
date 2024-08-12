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
using Interpolations

function edes_delay!(du, u, p, t, c_)
    # model equations with variable delay

    @unpack p_fixed, c, input = c_

    Mg, Gpl, Ipl, Gi = u
    
    mgmeal = p_fixed[10] * p[1]^p_fixed[10] * (t-input[3])^(p_fixed[10] - 1) * exp(-(p[1] * (t-input[3]))^p_fixed[10]) * input[1]
    mgpl = p_fixed[2] * Mg
    
    # glucose in gut
    du[1] = mgmeal - mgpl

    c11 = 0.043 * (p_fixed[11] + p_fixed[12]) / p_fixed[12]
    du[4] = (1.0/p[4]) * (Gpl - Gi)
    gliv = c[3] - p_fixed[3] * (Gpl - p_fixed[12]) - p_fixed[4] * c[4] * (Ipl - p_fixed[13])
    ggut = p_fixed[2] * (c[1] / (c[2] * input[2])) * Mg
    gnonit = c11 * (Gpl / (p_fixed[11] + Gpl)) 
    git = p_fixed[end] * c[4] * Ipl * (Gpl / (p_fixed[11] + Gpl)) # changed to p_fixed[6] as this reflect where the k6 is now
    gren = abs(((c[10] / (c[2] * input[2]) * (Gpl - c[8]))) * tanh(100 * (Gpl - c[8])) / 2 + (c[10] / (c[2] * input[2]) * (Gpl - c[8])) / 2)
    
    # glucose in plasma
    du[2] = gliv + ggut - gnonit - git - gren 

    ipnc = (c[4]^-1) * (p_fixed[end] * (Gpl - p_fixed[12]) + (p_fixed[7] / c[5]) * p_fixed[12] + (p_fixed[8] * c[6]) * du[2])
    iliv = c[12] * Ipl
    iif = p_fixed[9] * (Ipl - p_fixed[13])
    
    # insulin in plasma
    du[3] = ipnc - iliv - iif
    nothing
end

function edes_delayfix!(du, u, p, t, c_)
    # model equations
    @unpack p_fixed, c, input = c_

    Mg, Gpl, Ipl, Gi = u

    # Print initial state
    println("Time: $t")
    println("Initial state: Mg=$Mg, Gpl=$Gpl, Ipl=$Ipl, Gi=$Gi")
    println("Parameters: p=$p")
    println("Fixed parameters: p_fixed=$p_fixed")
    
    # Check for extreme values
    if any(x -> abs(x) > 1e10 || abs(x) < 1e-10, u)
        println("Warning: Extreme values detected in state variables")
    end

    mgmeal = p_fixed[10] * p[1]^p_fixed[10] * t^(p_fixed[10] - 1) * exp(-(p[1] * t)^p_fixed[10]) * input[1]
    mgpl = p_fixed[2] * Mg

    # glucose in gut
    du[1] = mgmeal - mgpl

    c11 = 0.043 * (p_fixed[11] + p_fixed[12]) / p_fixed[12] #- p[2] * p_fixed[13]
    du[4] = (1.0/2.5) * (Gpl - Gi)
    gliv = c[3] - p_fixed[3] * (Gpl - p_fixed[12]) - p_fixed[4] * c[4] * (Ipl - p_fixed[13])
    ggut = p_fixed[2] * (c[1] / (c[2] * input[2])) * Mg
    gnonit = c11 * (Gpl / (p_fixed[11] + Gpl)) 
    git = p_fixed[end] * c[4] * Ipl * (Gpl / (p_fixed[11] + Gpl)) # changed to p_fixed[6] 
    gren = abs(((c[10] / (c[2] * input[2]) * (Gpl - c[8]))) * tanh(100 * (Gpl - c[8])) / 2 + (c[10] / (c[2] * input[2]) * (Gpl - c[8])) / 2)
   
    # glucose in plasma    
    du[2] = gliv + ggut - gnonit - git - gren

    ipnc = (c[4]^-1) * (p_fixed[end] * (Gpl - p_fixed[12]) + (p_fixed[7] / c[5]) * p_fixed[12] + (p_fixed[8] * c[6]) * du[2])
    iliv = c[12] * Ipl
    iif = p_fixed[9] * (Ipl - p_fixed[13])

    # insulin in plasma
    du[3] = ipnc - iliv - iif

    # Print final derivatives
    println("Derivatives: dMg=$(du[1]), dGpl=$(du[2]), dIpl=$(du[3]), dGi=$(du[4])")
    
    # Check for extreme derivatives
    if any(x -> abs(x) > 1e10 || abs(x) < 1e-10, du)
        println("Warning: Extreme values detected in derivatives")
    end

    # Print intermediate calculated values
    println("Intermediate values:")
    println("mgmeal=$mgmeal, mgpl=$mgpl")
    println("gliv=$gliv, ggut=$ggut, gnonit=$gnonit, git=$git, gren=$gren")
    println("ipnc=$ipnc, iliv=$iliv, iif=$iif")

    # Check fixed k6 value
    println("Fixed k6 value: $(p_fixed[end])")
    if p_fixed[end] < 1e-10 || p_fixed[end] > 1e10
        println("Warning: k6 value may be out of expected range")
    end
    
    nothing
end

function loss(p, prob, data)
    data_glu = data[1,findall(!ismissing, data[1,4:end])]
    data_ins = data[2,findall(!ismissing, data[2,4:end])]

    ts_glu = parse.(Float64, names(data_glu))
    ts_ins = parse.(Float64, names(data_ins))

    prob = remake(prob)
    sol = solve(prob, RadauIIA5(), p=p, saveat=ts_glu; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
    
    loss = 0.0
    if any((s.retcode != :Success for s in sol))
        loss = Inf
    else
        loss_g = sum(abs2, ((sol[2,:] - Vector(data_glu))/maximum(Vector(data_glu))))
        loss_i = sum(abs2, ((sol[3,:] - Vector(data_ins))/maximum(Vector(data_ins))))
    
        loss = 100*(loss_i + loss_g)
    end
    return loss
end

function set_conditions(test, fixed_k6)
    if test == "OGTT"
        p_adj = [
            1.35e-2,  #k1
            3.80e-3,  #k5
        ]

        lb = [1e-7, 1e-7]
        ub = [0.3, 1.0]
    else
        p_adj = [
            1.35e-2,  #k1
            3.80e-3,  #k5
            2.5       #tau_g (delay)
        ]

        lb = [1e-7, 1e-7, 1e-10]
        ub = [0.3, 1.0, 20.0]
    end
    
    p_fixed = [
        1.35e-2,  #k1
        6.33e-1,  #k2
        5.00e-5,  #k3
        1.00e-3,  #k4
        3.80e-3,  #k5
        1.15,     #k7
        4.71,     #k8
        1.08e-2,  #k9
        1.35,     #sigma /index 10
        0.63,     #Km /index 11
        5.0,      #Gb /index 12
        8.0,      #Ib /index 13
        5.0,      #gpllowerbound
        2.84e-1,
        fixed_k6 #k6 (now fixed)
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
    
    tspan = (0.0, 240.0)    
    u0 = [0.0, p_fixed[12], p_fixed[13], p_fixed[12]]

    return p_adj, lb, ub, p_fix, u0, tspan
end

function update_conditions!(p_adj, p_fix, u0, data; test="OGTT")
    # used to update the conditions in set_condition to the particular data used in calibration

    u0[2] = data[1, 4] # Glucose, tp=0
    u0[3] = data[2, 4] # Insulin, tp=0
    u0[4] = data[1, 4] # interstitial glucose tp=0

    p_fix.p_fixed[12] = data[1, 4]  # update to measured basal glucose
    p_fix.p_fixed[13] = data[2, 4]  # update to measured basal insulin

    return u0
end

function PLA(opt, lb, ub, prob, data, storage; test)
    if test == "OGTT"
        scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100)]
        theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200)]
    else
        scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100), (1e-9, 25.0)]
        theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200), (1e-10, 30.0)]
    end
    
    α = loss(opt, prob, data) + 3.84

    num_params = length(opt);

    intervals_ = Vector{ParamInterval}(undef, num_params);

    for i in 1:num_params
        @time intervals_[i] = get_interval(
            opt,
            i,
            x->loss(x,prob,data),
            :CICO_ONE_PASS,
            loss_crit= α,
            theta_bounds=theta_bounds,
            scan_bounds=scan_bounds[i],
            scale=fill(:direct, length(opt)),
            scan_tol=1e-5,
            local_alg = :LN_NELDERMEAD,
            max_iter=10^5
        );
        storage[1, i*3-2] = intervals_[i].result[1].value 
        storage[1, i*3] = opt[i]  # Store the optimized value
        storage[1, i*3+1] = intervals_[i].result[2].value
    end;

    # Store the fixed k6 value
    if test == "OGTT"
        storage[1, 7] = prob.p[end]
    else
        storage[1, 10] = prob.p[end]
    end

    return
end

function bands(opt, sol, lb, ub, prob, data, ts_sim; test="OGTT")
    function bands_g(p, i)
        prob = remake(prob) 
        sol_ = solve(prob, RadauIIA5(), p=p, saveat=ts_sim; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        return sol_[2,i]
    end

    function bands_sg(p, i)
        prob = remake(prob) 
        sol_ = solve(prob, RadauIIA5(), p=p, saveat=ts_sim; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        return sol_[4,i]
    end

    function bands_i(p, i)
        prob = remake(prob) 
        sol_ = solve(prob, RadauIIA5(), p=p, saveat=ts_sim; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
        return sol_[3,i]
    end
    
    if test == "OGTT"
        scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100)]
        theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200)]
    else
        scan_bounds = [(1e-10, ub[1]*100), (1e-10, ub[2]*100), (1e-9, 25.0)]
        theta_bounds = [(0.0, ub[1]*200), (0.0, ub[2]*200), (1e-10, 30.0)]
    end

    α = loss(opt, prob, data) + 3.84

    intervals_g = Vector{ParamInterval}(undef,length(ts_sim))
    if test == "CGM"
        for i in eachindex(ts_sim)
            @time intervals_g[i] = get_interval(
                opt,
                p->bands_sg(p,i),
                x->loss(x,prob,data),
                :CICO_ONE_PASS,
                loss_crit = α,
                theta_bounds=theta_bounds,
                scan_bounds=(0.00001,500.0), 
                scan_tol = 1e-5,
                local_alg = :LN_NELDERMEAD,
            )
        end
    else
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
    end

    glbands = [iv.result[1].value for iv in intervals_g]
    gubands = [iv.result[2].value for iv in intervals_g]

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

    ilbands = [iv.result[1].value for iv in intervals_i]
    iubands = [iv.result[2].value for iv in intervals_i]
    return glbands, gubands, ilbands, iubands
end

function plot_bands(sol, data, glbands, gubands, ilbands, iubands; test)
    if test == "CGM"
        u_g = sol[4]
    else
        u_g = sol[2]
    end
    u_i = sol[3]

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
end

function gen_results(;test="OGTT")
    println("Starting gen_results function")
    data = CSV.read("test.csv", DataFrame)
    
    println("Data structure:")
    println(first(data, 5))
    
    sort!(data, :time)
    data_wide = unstack(data, [:ID, :metab, :Condition], :time, :VAL)
    
    println("\nWide data structure:")
    println(first(data_wide, 5))

    grp_data = groupby(data_wide, [:ID, :Condition])
    println("\nNumber of groups: ", length(grp_data))

    # Calculate average k6
    average_k6 = 1.167872529    # This is the average for k6_test. Replace with actual average calculation if needed.
    println("Average k6 value: ", average_k6)

    simulated = copy(data_wide[:,1:3])
    sim_ts = DataFrame(missings(Float64,nrow(simulated),length(names(data_wide))-3), :auto)
    simulated = hcat(simulated, sim_ts)
    glu_bands = copy(simulated)
    ins_bands = copy(simulated)
    grp_simulated = groupby(simulated, [:ID, :Condition])
    grp_glu_bands = groupby(glu_bands, [:ID, :Condition])
    grp_ins_bands = groupby(ins_bands, [:ID, :Condition])

    estimates = copy(unique(data_wide[:,[:ID, :Condition]]))
    if test == "OGTT"
        est_p = DataFrame(Array{Union{Nothing,Float64}}(nothing,nrow(estimates),6), [:k1_l,:k1,:k1_r, :k5_l,:k5,:k5_r])
    else
        est_p = DataFrame(Array{Union{Nothing,Float64}}(nothing,nrow(estimates),9), [:k1_l,:k1,:k1_r, :k5_l,:k5,:k5_r, :tauG_l,:tauG,:tauG_r])
    end

    estimates = hcat(estimates, est_p)
    estimates.fixed_k6 = fill(average_k6, nrow(estimates))  # Add fixed k6 column
    grp_estimates = groupby(estimates, [:ID, :Condition])

    residuals = copy(data_wide)
    residuals[:,(4:end)] .= missing
    grp_residuals = groupby(residuals, [:ID, :Condition])

    p_adj, lb, ub, p_fix, u0, tspan = set_conditions(test, average_k6)
    
    n = length(grp_data)

    for i in 1:n
        try
            println("\nProcessing group $i of $n")
            println(grp_data[i][1,:ID], grp_data[i][1,:Condition])

            u0 = update_conditions!(p_adj, p_fix, u0, grp_data[i]; test)
            if test == "OGTT"
                prob = ODEProblem((du, u, p, t) -> edes_delayfix!(du, u, p, t, p_fix), u0, tspan, p_adj)
            else
                prob = ODEProblem((du, u, p, t) -> edes_delay!(du, u, p, t, p_fix), u0, tspan, p_adj)
            end

            # Check if ODE can be solved with initial parameters
            try
                sol = solve(prob, RadauIIA5(), abstol=1e-8, reltol=1e-6)
                if sol.retcode != :Success
                    println("ODE solver failed with initial parameters. Retcode: ", sol.retcode)
                    continue
                end
            catch e
                println("Error in initial ODE solve:")
                println(e)
                continue
            end

            # Check if loss function is working
            try
                initial_loss = loss(p_adj, prob, grp_data[i])
                println("Initial loss: ", initial_loss)
            catch e
                println("Error in initial loss calculation:")
                println(e)
                continue
            end

            optf = OptimizationFunction((x, p_adj) -> loss(x, prob, grp_data[i])) 
            opt_prob = OptimizationProblem(optf, p_adj, lb=lb, ub=ub)
            
            println("Initial conditions: ", u0)
            println("Lower bounds: ", lb)
            println("Upper bounds: ", ub)
            println("Starting optimization...")
            
            opt_sol = solve(opt_prob, MultistartOptimization.TikTak(1000), NLopt.LN_NELDERMEAD(), abstol = 1e-7, reltol = 1e-7) 
            println("Optimization complete. Result: ", opt_sol)

            _prob = remake(prob, p=opt_sol.u)
            sol_opt = solve(_prob, RadauIIA5(), saveat=0:240; abstol = 1e-7, reltol = 1e-5, maxiters=1e7)
            sol = Array(sol_opt)

            ts_data = parse.(Int, names(grp_data[i])[4:end])
            println("Time points: ", ts_data)

            interp_sol = [LinearInterpolation(0:240, sol[j,:]) for j in 1:4]

            println("Dimensions of sol: ", size(sol))
            println("Length of ts_data: ", length(ts_data))

            if test == "OGTT"
                grp_simulated[i][1,4:end] = interp_sol[2].(ts_data)
                grp_estimates[i][1,[:k1, :k5]] = opt_sol.u
            else
                grp_simulated[i][1,4:end] = interp_sol[4].(ts_data)
                grp_estimates[i][1,[:k1, :k5, :tauG]] = opt_sol.u
            end
            grp_simulated[i][2,4:end] = interp_sol[3].(ts_data)

            println("Dimensions of grp_simulated[i]: ", size(grp_simulated[i]))

            gluc_data = Array(grp_data[i][isequal.(grp_data[i].metab,"gluc"), 4:end])[1,:]
            ins_data = Array(grp_data[i][isequal.(grp_data[i].metab,"ins"), 4:end])[1,:]
            
            println("Dimensions of gluc_data: ", size(gluc_data))
            println("Dimensions of ins_data: ", size(ins_data))

            if test == "OGTT"
                grp_residuals[i][1, 4:end] = gluc_data - Vector(grp_simulated[i][1,4:end])
            else
                grp_residuals[i][1, 4:end] = gluc_data - Vector(grp_simulated[i][1,4:end])
            end
            grp_residuals[i][2, 4:end] = ins_data - Vector(grp_simulated[i][2,4:end])

            println("Dimensions of grp_residuals[i]: ", size(grp_residuals[i]))

            glbands, gubands, ilbands, iubands = bands(opt_sol.u, sol_opt, lb, ub, _prob, grp_data[i], Float64.(ts_data); test)
            grp_glu_bands[i][1,4:end] = glbands
            grp_glu_bands[i][2,4:end] = gubands
            grp_ins_bands[i][1,4:end] = ilbands
            grp_ins_bands[i][2,4:end] = iubands
            
            try
                plot_bands(interp_sol, grp_data[i], glbands, gubands, ilbands, iubands; test)
            catch e
                println("Error in plot_bands: ", e)
            end
        catch e
            println("Error processing group $i:")
            println(e)
            println(stacktrace(catch_backtrace()))
        end
    end

    println("Processing complete. Cleaning up data...")

    dropmissing!(simulated, :x1)
    dropmissing!(glu_bands, :x1)
    dropmissing!(ins_bands, :x1)
    
    println("Returning results...")

    return simulated, estimates, residuals, glu_bands, ins_bands
end