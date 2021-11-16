# Weighted estimators of survival under a generalized representative intervention for organ transplant

This repository contains the code to simulate survival data for patients awaiting organ transplant and fit a time-varying weighted Kaplan-Meier estimator of survival under a generalized representative intervention, as introduced by Lyden, et al. in "Transportability of causal inference under random dynamic treatment regimes for kidney-pancreas transplantation."

A generalized representative intervention is a special case of a random dynamic treatment regime that **assigns treatment by drawing from the distribution of treatments that would have been received by identical patients in a target population, if those patients had been intervened on to comply with the rules of the regime**. 

For the simulation study, the strategy of interest is "wait for deceased-donor transplant without receiving living-donor transplant," and the target population is a one of two transplant centers in the study population. Both the rate of transplant and the quality of offered organs are more favorable at the target center.

The getRobustCI.R script includes a function to compute a sandwich variance estimator for log survival, which is not currently available from the survival package in R, and a second function to construct 95% pointwise confidence intervals using the log transformation.
