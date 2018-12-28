# task_complexity_2018
Simulation code for the computational model of arbitration and fMRI results

[1] Simulation Code : main.m
This will generate a file for each subject in modelRLsource/result_simul. Every simulation file(.mat) contains SBJ structure variable. SBJ{1, 1}.model_BayesArb.param is the optimized free parameters, and SBJ{1, 1}.model_BayesArb.val is sum of negative log likelihood for each subject. 

The sum of negative log likelihood further used to calculate BIC and model comparison (RFX BMS).


[2] Group effect fMRI results : ./group_effect_fmri
