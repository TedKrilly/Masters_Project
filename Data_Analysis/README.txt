Glucose and Insulin Data Analysis Scripts

R scripts for analysing glucose data from clinical trials. Handles CGM, OGTT, and HFMM data processing
and analysis. The analysis portion of my masters project working with clinical trial data.

The Content: 
1. Masters_Project.Rmd : Main script comparing OGTT/CGM fits with density plots. Almost all the plots 
and analysis are in this file. 

2. Baseline_Model_Data_Processing.Rmd: Cleans up CGM & HFMM data

3. CGM_basic_metrics.Rmd: Basic CGM stats and correlations

4. HFMM_Quality_Of_Fit.Rmd: Checks how well the models fit -> Omitted from final project

5. k5_analysis.Rmd: Compares k5 values between different tests -> Comparison of the different k5 values

Issues : 
1. File picker is janky, you'll need to select files manually
2. Dates are a mess, formats vary between datasets
3. Some missing data - scripts handle it but check your outputs
4. k5 analysis assumes specific column names, check before running
