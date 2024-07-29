using Pumas

edes_model = @emmodel begin
    @param begin
        k1 ~ 1 | LogNormal
        k5 ~ 1 | LogNormal
        k6 ~ 1 | LogNormal
        τ_g ~ 1 | LogNormal
    end

    @random begin
        η_k1 ~ 1 | Normal
        η_k5 ~ 1 | Normal
        η_k6 ~ 1 | Normal
        η_τg ~ 1 | Normal
    end

    @covariates begin
        weight
        data_type  # 'OGTT' or 'CGM'
    end

    @pre begin
        # Fixed parameters
        p_fixed = [1.35e-2, 6.33e-1, 5.00e-5, 1.00e-3, 3.80e-3,
                   5.82e-1, 1.15, 4.71, 1.08e-2, 1.35,
                   0.63, 5.0, 8.0, 5.0, 2.84e-1]
        
        # Constants
        c = [0.005551, 17/70, 0.043, 1.0, 31.0,
             3.0, 13.0/70.0, 9.0, 30.0, 0.1,
             0.043 * (5.63) / 5.0 - 3.80e-3 * 1 * 8.0,
             1.15 * (5.0 / (31 * 8.0))]
        
        # Input parameters (these should be provided as covariates or in the event data)
        input = [75000.0, weight, 0.0]  # meal size, body weight, meal start time

        # Derived constants
        c11 = 0.043 * (p_fixed[11] + p_fixed[12]) / p_fixed[12]

        # Apply random effects
        k1_ind = k1 * exp(η_k1)
        k5_ind = k5 * exp(η_k5)
        k6_ind = k6 * exp(η_k6)
        τ_g_ind = τ_g * exp(η_τg)
    end

    @init begin
        Mg = 0.0
        Gpl = p_fixed[12]  # Basal glucose level
        Ipl = p_fixed[13]  # Basal insulin level
        Gi = p_fixed[12]   # Initial interstitial glucose
    end
    
    @dynamics begin
        mgmeal = p_fixed[10] * k1_ind^p_fixed[10] * (t-input[3])^(p_fixed[10] - 1) * exp(-(k1_ind * (t-input[3]))^p_fixed[10]) * input[1]
        mgpl = p_fixed[2] * Mg

        Mg' = mgmeal - mgpl

        # Always calculate Gi for CGM, but use it only when data_type is CGM
        Gi' = (1.0/τ_g_ind) * (Gpl - Gi)

        gliv = c[3] - p_fixed[3] * (Gpl - p_fixed[12]) - p_fixed[4] * c[4] * (Ipl - p_fixed[13])
        ggut = p_fixed[2] * (c[1] / (c[2] * input[2])) * Mg
        gnonit = c11 * (Gpl / (p_fixed[11] + Gpl))
        git = k5_ind * c[4] * Ipl * (Gpl / (p_fixed[11] + Gpl))
        gren = abs(((c[10] / (c[2] * input[2]) * (Gpl - c[8]))) * tanh(100 * (Gpl - c[8])) / 2 + (c[10] / (c[2] * input[2]) * (Gpl - c[8])) / 2)

        Gpl' = gliv + ggut - gnonit - git - gren

        ipnc = (c[4]^-1) * (k6_ind * (Gpl - p_fixed[12]) + (p_fixed[7] / c[5]) * p_fixed[12] + (p_fixed[8] * c[6]) * Gpl')
        iliv = c[12] * Ipl
        iif = p_fixed[9] * (Ipl - p_fixed[13])

        Ipl' = ipnc - iliv - iif
    end

    @derived begin
        glucose := data_type == "CGM" ? Gi : Gpl
    end

    @error begin
        glucose ~ data_type == "CGM" ? Normal(glucose, 0.1 * glucose) : ProportionalNormal(glucose)
    end
end