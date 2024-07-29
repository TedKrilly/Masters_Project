Pumas In Built Functions

Pumas uses in built optimisation algorithms to estimate the model parameters by minimising the difference
between the observed and predicted data
  As a result dont need to include Edes model loss function - dont know how Pumas optimiser works

Omitted Plotting functions and PLA as for the moment I am interested in producing parameter estimates exclusively

Model can handle both OGTT and CGM data by adjusting the derived variables and error model based on 
the data_type covariate

Assumptions of the model:

    The model assumes that the parameters, random effects, and input parameters are log-normally distributed.
    The model assumes that the error in the glucose measurements is normally distributed when data_type is "CGM" and proportional to the glucose measurements when data_type is not "CGM".

Pumas-Func-Logic

@emmodel macro is used as SAEM estimation requires this model form

@param block: used to define model parameters
    - four key parameters (k1, k5, k6, τ_g) 
    - assume LogNormal distributions

@random block: 

@covariates:

@preprocessing: pre block is used for preprocessing before the model simulation.
p_fixed contains fixed parameter values.
c contains constants used in the model.
input contains input parameters like meal size, body weight, and meal start time.
c11 is a derived constant.
k1_ind, k5_ind, k6_ind, and τ_g_ind are individual-specific parameters adjusted by random effects.

@init block: defines the initial conditions for the state variables
Mg, Gpl, Ipl, and Gi are the initial values for meal glucose, plasma glucose, plasma insulin, 
and interstitial glucose, respectively.

@dynamics: define the EDES model differential equations that determine the model dynamics

@derived block: 
glucose is set to Gi if data_type is "CGM", otherwise it is set to Gpl.

@error block: If data_type is "CGM", the glucose measurement error is modeled as a normal distribution with a standard deviation proportional to the glucose level.
Otherwise, a proportional normal error model is used.
