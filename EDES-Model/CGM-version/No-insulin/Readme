This EDES-Model variant was created to investigate what happens when you try to estimate k1,k5, and k6 when you have only CGM glucose (no insulin). 
I would expect that we can't find unique solutions( this is something called non-identfiabiliy) 

By making these changes, we're attempting to estimate k1, k5, and k6 using only CGM glucose data. This is likely to lead to identifiability issues because:

k1 influences how quickly glucose appears in the blood after a meal.
k5 relates to insulin-dependent glucose uptake.
k6 relates to insulin secretion.

Without insulin data, it becomes difficult to distinguish between the effects of these parameters. For example:

A high glucose peak could be due to rapid glucose appearance (high k1) or poor insulin-dependent uptake (low k5).
A rapid decline in glucose could be due to high insulin sensitivity (high k5) or rapid insulin secretion (high k6).


Expected outcomes:

When running this modified version, we expect to see:

The optimization algorithm may struggle to converge on unique solutions.
Parameter estimates may hit the bounds of their allowed ranges.
Confidence intervals for parameter estimates may be very wide.
Different initial guesses may lead to wildly different final estimates.
The model may fit the glucose data reasonably well, but with unrealistic parameter values.

These outcomes would demonstrate that the model parameters are not identifiable using CGM glucose data alone, highlighting the importance of including insulin measurements in the parameter estimation process.
This analysis serves to illustrate the limitations of using glucose data alone and the value of including insulin measurements in model calibration for accurate estimation of insulin sensitivity and secretion parameters.
