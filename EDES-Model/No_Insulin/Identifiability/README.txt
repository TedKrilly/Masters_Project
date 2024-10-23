Parameter Identifiability Analysis Variant
This version of the streamlined glucose model adds comprehensive parameter identifiability
analysis. It focuses specifically on CGM data and provides detailed assessment of
parameter identifiability through profile likelihood analysis (PLA).
Key Features:

Enhanced PLA Function


Tracks identifiability status for each parameter
Categorizes parameters as:

Identifiable
Practically unidentifiable
Structurally unidentifiable


Saves detailed results to CSV file
Includes error tracking and reporting
Parameters analyzed: k1, k5, k6, tauG


Added Outputs


pla_results_CGM_noins.csv: Detailed PLA results including

Parameter point estimates
Confidence intervals
Identifiability classification
Error information


CGMestimates_noins.csv: Parameter estimates
CGMsimulates_noins.csv: Model simulations
CGMresiduals_noins.csv: Fit residuals


Identifiability Classification


Based on confidence interval width relative to parameter range
< 10% of parameter range: Identifiable
< 100% of parameter range: Practically unidentifiable


100% of parameter range: Structurally unidentifiable




Error Handling


Tracks failed PLA calculations
Stores error messages for debugging
Continues analysis even if individual parameters fail
Reports comprehensive error summaries
