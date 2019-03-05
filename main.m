% 1. model fitting 
% set the job_opt.option_optimizing_model as 0 in the
% 'batch_model_cog_regressor_gen_v5_indi_complex' line 49
for sub = 1:1:24
    batch_model_cog_regressor_gen_v5_indi_complex(sub, 0, 0, 10)
end

% 2. Generating regressors
% set the job_opt.option_optimizing_model as 1 in the
% 'batch_model_cog_regressor_gen_v5_indi_complex' line 49
batch_model_cog_regressor_gen_v5_indi_complex(sub, 0, 0, 10)
