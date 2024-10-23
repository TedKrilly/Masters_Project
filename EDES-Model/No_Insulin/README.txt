This EDES-Model variant was created to investigate what happens when you try to estimate k1,k5, and k6 when you have only CGM glucose (no insulin). 
I would expect that we can't find unique solutions( this is something called non-identfiabiliy) 

By making these changes, we're attempting to estimate k1, k5, and k6 using only CGM glucose data. This is likely to lead to identifiability issues because:

k1 influences how quickly glucose appears in the blood after a meal.
k5 relates to insulin-dependent glucose uptake.
k6 relates to insulin secretion.

Without insulin data, it becomes difficult to distinguish between the effects of these parameters. For example:

A high glucose peak could be due to rapid glucose appearance (high k1) or poor insulin-dependent uptake (low k5).
A rapid decline in glucose could be due to high insulin sensitivity (high k5) or rapid insulin secretion (high k6).


This is a modified version of the glucose-insulin model that focuses solely on glucose dynamics.
The key modification is the removal of insulin fitting from the loss function, making the
model more robust when insulin data is incomplete or unreliable.
Key Modifications from Original Model:

Loss Function


Original weighted both glucose and insulin fits
Modified to only consider glucose error: loss = 100*loss_g
Still tracks insulin dynamics but doesn't fit to insulin data
More stable optimization due to simpler objective function


Data Requirements


Primary input: glucose measurements
Insulin measurements optional
Works with both OGTT and CGM data
Added robustness checks for data completeness


Error Handling


Checks for missing data before processing
Verifies minimum data points for CGM analysis
Gracefully skips problematic datasets
Clear console messages for skipped analyses

Input Data Format (test.csv):
Required columns:

ID: Subject identifier
metab: Measurement type (gluc/ins)
Condition: Test condition
test: OGTT or CGM
time: Measurement timestamps
VAL: Measured values

Output Files:

OGTTestimates.csv: Parameter estimates for OGTT data
OGTTsimulates.csv: Model simulations
OGTTresiduals.csv: Fit residuals
Plus CGM equivalents if CGM data present
Individual plots for each subject/condition
