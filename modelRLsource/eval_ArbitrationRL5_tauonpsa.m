function [data_out] = eval_ArbitrationRL5_tauonpsa(param_in, data_in, mode)
% Load the parameter files from the optimization folder, and then evaluate the model.



%% create MAP
map_opt.transition_prob_seed=[0.5 0.5];
% map_opt.reward_seed=[40 20 10 0];
% [myMap N_state N_action N_transition]=Model_Map_Init2('sangwan2012b',map_opt);
% [myMap_new N_state N_action N_transition]=Model_Map_Init2('sangwan2012c',map_opt);

map_opt.transition_prob_seed=[0.9 0.1];
[myMap_new N_state N_action N_transition]=Model_Map_Init2('sangwan2014b',map_opt);
myMap=myMap_new;


%% create my arbitrator
myArbitrator=Bayesian_Arbitration_Init_v2(N_state,N_action,N_transition,mode);


%% create my RL
myState=Model_RL_Init_v2(N_state,N_action,N_transition,mode);



%% model parameter - for the functional mode of RL
% SARSA model
param_sarsa.gamma=1.0; % fixed - not actual parameter
%     param_sarsa.alpha=0.15; % learning rate (0.1~0.2)
%     param_sarsa.tau=0.5; % decision focus parameter
% FWD model (for the latent stage)
% param_fwd0.alpha=0.1; % learning rate
%     param_fwd0.tau=param_sarsa.tau; % decision focus parameter
% FWD model
%     param_fwd.alpha=0.15; % learning rate
%     param_fwd.tau=param_sarsa.tau; % decision focus parameter





%% SIMULATION


%% parameter plug-in
% param_in(1): myArbitrator.PE_tolerance_m1 (m1's threshold for zero PE)
% param_in(2): myArbitrator.PE_tolerance_m2 (m2's threshold for zero PE)
% param_in(3): myArbitrator.A_12
% param_in(x): myArbitrator.B_12 : based on A12
% param_in(4): myArbitrator.A_21
% param_in(x): myArbitrator.B_21 : based on A21
% param_in(5): myArbitrator.tau_softmax/param_sarsa.tau/param_fwd.tau : better to fix at 0,2. This should be determined in a way that maintains softmax values in a reasonable scale. Otherwise, this will drive the fitness value!
% param_in(6): % param_sarsa.alpha/param_fwd.alpha 0.01~0.2 to ensure a good "state_fwd.T" in phase 1
param_fixed(1)=1; % 1: fwd-start, 2:sarsa-start
param_fixed(2)=1; % myArbitrator.p
param_fixed(3)=1e-1; % myArbitrator.Time_Step : time constant (1e0 (fast) ~ 1e-2 (slow))

pop_id=1; % fixed
Sum_NegLogLik=0.0;
% arbitrator
myArbitrator.PE_tolerance_m1=param_in(pop_id,1); % defines threshold for zero PE
% myArbitrator.PE_tolerance_m2=param_in(pop_id,2); % defines threshold for zero PE
myArbitrator.m2_absPEestimate_lr=param_in(pop_id,2); % defines the learning rate of RPE estimator
myArbitrator.Time_Step=param_fixed(3); % the smaller, the slower
switch mode.param_length
    case 6
        myArbitrator.A_12=param_in(pop_id,3);
        myArbitrator.B_12=log(myArbitrator.A_12/mode.boundary_12-1);
        myArbitrator.A_21=param_in(pop_id,4);
        myArbitrator.B_21=log(myArbitrator.A_21/mode.boundary_21-1);
        % SARSA
        param_sarsa.alpha=param_in(pop_id,6); % learning rate (0.1~0.2)
        % FWD
        param_fwd.alpha=param_in(pop_id,6);
        % Softmax parameter for all models
        myArbitrator.tau_softmax=param_in(pop_id,5); % use the same value as sarsa/fwd
        param_sarsa.tau=param_in(pop_id,5);
        param_fwd.tau=param_in(pop_id,5);
        %     case 8 % ocg model
        %         myArbitrator.A_12=param_in(pop_id,3);
        %         myArbitrator.B_12=param_in(pop_id,4);
        %         myArbitrator.A_21=param_in(pop_id,5);
        %         myArbitrator.B_21=param_in(pop_id,6);
        %         % SARSA
        %         param_sarsa.alpha=param_in(pop_id,8); % learning rate (0.1~0.2)
        %         % FWD
        %         param_fwd.alpha=param_in(pop_id,8);
        %         % Softmax parameter for all models
        %         myArbitrator.tau_softmax=param_in(pop_id,7); % use the same value as sarsa/fwd
        %         param_sarsa.tau=param_in(pop_id,7);
        %         param_fwd.tau=param_in(pop_id,7);
    case 8 % cog model
        myArbitrator.A_12=param_in(pop_id,3);
        myArbitrator.B_12=log(myArbitrator.A_12/mode.boundary_12-1);
        myArbitrator.C_12=param_in(pop_id,4);
        myArbitrator.A_21=param_in(pop_id,5);
        myArbitrator.B_21=log(myArbitrator.A_21/mode.boundary_21-1);
        myArbitrator.C_21=param_in(pop_id,6);
        % SARSA
        param_sarsa.alpha=param_in(pop_id,8); % learning rate (0.1~0.2)
        % FWD
        param_fwd.alpha=param_in(pop_id,8);
        % Softmax parameter for all models
        myArbitrator.tau_softmax=param_in(pop_id,7); % use the same value as sarsa/fwd
        param_sarsa.tau=param_in(pop_id,7);
        param_fwd.tau=param_in(pop_id,7);
end

% arbitrator start
myArbitrator.ind_active_model=param_fixed(1);
if(myArbitrator.ind_active_model==1)
    myArbitrator.m1_prob_prev=0.7; % do not use 0.5 which causes oscillation.
    myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
else
    myArbitrator.m1_prob_prev=0.3;
    myArbitrator.m2_prob_prev=1-myArbitrator.m1_prob_prev;
end
% non-linear weight : p
myArbitrator.p=param_fixed(2);

% arbitrator functioning mode
myArbitrator.opt_ArbModel=mode.opt_ArbModel;
myArbitrator.opt_cogload=mode.opt_cogload;

%% Devaluation/coningency degradation/risk schedule
% % (1) devaluation schedule
% pt_devaluation=[1000];%[41:1:80];%[1:1:40]; % second phase 1~80. if devaluation point > num_max_trial, then no devaluation
% pt_ref=pt_devaluation(1);
% % (2) degradation schedule
% pt_degradation=[1000];
% % (3) transition prob schedule
% pt_prob_chg_max_risk=[1000];%[1:1:30 61:1:90 121:1:150];%[11:1:40]; % second phase 1~80.
% pt_prob_chg_min_risk=[1000];%[31:1:60 91:1:120 151:1:180];%[41:1:80]; %[41:1:80]; % second phase 1~80.
% pt_prob_chg_reversed=[1000]; % second phase 1~80.
% % (4) reward prob schedule
% pt_prob_rwd_chg=[1000];
% % (5) goal-directed mode (keep changing rwd values)
% pt_devaluation_goal=[1000];%[31:1:60 91:1:120 151:1:180];
% rwd_aversion_factor=0; % 1: no change, 0: zero reward (ok), -1:opposite (best)
%
% % scenario 3
% ss0=[1:1:30]; ss1=[31:1:60]; ss2=[61:1:90]; ss3=[91:1:120];   ss4=[121:1:150];   ss5=[151:1:180];
% pt_prob_chg_max_risk=[ss0 ss2 ss4];
% pt_prob_chg_min_risk=[ss1 ss3 ss5];
% pt_devaluation_goal=[ss1 ss3 ss5];

pt_devaluation=1000;


for ind_dev=[1:1:length(pt_devaluation)]
    
    %% Multiple simulations
    win_fwd=0; win_sarsa=0;
    
    tot_num_sbj=size(data_in,2);
    for ll=1:1:tot_num_sbj % each subject
        
        
        n_tmp=0;
        if(mode.out==99) % in regressor generation mode, do simulation only once
            mode.total_simul=1;
        end
        
        
        for kk=1:1:mode.total_simul
            
            % regressor part : DO NOT change the order ############################
            
            if(mode.out==99)
                %                 disp(sprintf('- generating regressor for SBJ%02d/%02d...',ll,tot_num_sbj));
                data_in{1,ll}.regressor{1,1}.name='SPE';     data_in{1,ll}.regressor{1,1}.value=[]; % SPE
                data_in{1,ll}.regressor{1,1}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: SPE'};
                data_in{1,ll}.regressor{1,2}.name='RPE';     data_in{1,ll}.regressor{1,2}.value=[]; % RPE
                data_in{1,ll}.regressor{1,2}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: RPE'};
                data_in{1,ll}.regressor{1,3}.name='uncertaintyM1';     data_in{1,ll}.regressor{1,3}.value=[];
                data_in{1,ll}.regressor{1,3}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: uncertainty of fwd model'};
                data_in{1,ll}.regressor{1,4}.name='uncertaintyM2';     data_in{1,ll}.regressor{1,4}.value=[];
                data_in{1,ll}.regressor{1,4}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: uncertainty of sarsa model'};
                data_in{1,ll}.regressor{1,5}.name='meanM1';     data_in{1,ll}.regressor{1,5}.value=[];
                data_in{1,ll}.regressor{1,5}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: mean belief of fwd model'};
                data_in{1,ll}.regressor{1,6}.name='meanM2';     data_in{1,ll}.regressor{1,6}.value=[];
                data_in{1,ll}.regressor{1,6}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: mean belief of sarsa model'};
                data_in{1,ll}.regressor{1,7}.name='invFanoM1';     data_in{1,ll}.regressor{1,7}.value=[];
                data_in{1,ll}.regressor{1,7}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: invFano of fwd model'};
                data_in{1,ll}.regressor{1,8}.name='invFanoM2';     data_in{1,ll}.regressor{1,8}.value=[];
                data_in{1,ll}.regressor{1,8}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7-9: invFano of sarsa model'};
                data_in{1,ll}.regressor{1,9}.name='weigtM1';     data_in{1,ll}.regressor{1,9}.value=[];
                data_in{1,ll}.regressor{1,9}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: choice prob of fwd model'};
                data_in{1,ll}.regressor{1,10}.name='weigtM2';     data_in{1,ll}.regressor{1,10}.value=[]; % (1-WgtModel1)
                data_in{1,ll}.regressor{1,10}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: choice prob of sarsa model'};
                data_in{1,ll}.regressor{1,11}.name='Qfwd';     data_in{1,ll}.regressor{1,11}.value=[]; %
                data_in{1,ll}.regressor{1,11}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: Q(s,a) of fwd model'};
                data_in{1,ll}.regressor{1,12}.name='Qsarsa';     data_in{1,ll}.regressor{1,12}.value=[]; %
                data_in{1,ll}.regressor{1,12}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: Q(s,a) of sarsa model'};
                data_in{1,ll}.regressor{1,13}.name='Qarb';     data_in{1,ll}.regressor{1,13}.value=[]; %
                data_in{1,ll}.regressor{1,13}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: Q(s,a) of the arbitrator'};
                data_in{1,ll}.regressor{1,14}.name='dQbwdEnergy';     data_in{1,ll}.regressor{1,14}.value=[]; %
                data_in{1,ll}.regressor{1,14}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: sum of (Q(s,a) change).^2 in backward update of fwd model'};
                data_in{1,ll}.regressor{1,15}.name='dQbwdMean';     data_in{1,ll}.regressor{1,15}.value=[]; %
                data_in{1,ll}.regressor{1,15}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: norm of Q-vector change in backward update of fwd model'};
                data_in{1,ll}.regressor{1,16}.name='duncertaintyM1';     data_in{1,ll}.regressor{1,16}.value=[]; %
                data_in{1,ll}.regressor{1,16}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: change of M1 uncertainty'};
                data_in{1,ll}.regressor{1,17}.name='dinvFanoM1';     data_in{1,ll}.regressor{1,17}.value=[]; %
                data_in{1,ll}.regressor{1,17}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: change of M1 invFano'};
                data_in{1,ll}.regressor{1,18}.name='TR_alpha';     data_in{1,ll}.regressor{1,18}.value=[]; %
                data_in{1,ll}.regressor{1,18}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: transition rate for MF2MB'};
                data_in{1,ll}.regressor{1,19}.name='TR_beta';     data_in{1,ll}.regressor{1,19}.value=[]; %
                data_in{1,ll}.regressor{1,19}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: transition rate for MB2MF'};
                data_in{1,ll}.regressor{1,20}.name='dinvFano12';     data_in{1,ll}.regressor{1,20}.value=[]; % invF1-invF2 >> will be overwritten at the end for mean correction!
                data_in{1,ll}.regressor{1,20}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: invFanoM1 minus invFanoM2'};
                data_in{1,ll}.regressor{1,21}.name='ABSdinvFano12';     data_in{1,ll}.regressor{1,21}.value=[]; % |invF1-invF2|  >> will be overwritten at the end for mean correction!
                data_in{1,ll}.regressor{1,21}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: abs(invFanoM1 minus invFanoM2)'};
                data_in{1,ll}.regressor{1,22}.name='MAXinvFano12';     data_in{1,ll}.regressor{1,22}.value=[]; % max(invF1,invF2)  >> will be overwritten at the end for mean correction!
                data_in{1,ll}.regressor{1,22}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: max(invFanoM1,invFanoM2)'};
                data_in{1,ll}.regressor{1,23}.name='CONFLICTinvFano12';     data_in{1,ll}.regressor{1,23}.value=[]; % conflict btw invF1&invF2 >> will be overwritten at the end for mean correction!
                data_in{1,ll}.regressor{1,23}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: max(invFanoM1,invFanoM2)'};
                data_in{1,ll}.regressor{1,24}.name='Qchan';     data_in{1,ll}.regressor{1,24}.value=[]; % SPE
                data_in{1,ll}.regressor{1,24}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: Qchan'};
                data_in{1,ll}.regressor{1,25}.name='complexity';     data_in{1,ll}.regressor{1,25}.value=[]; % SPE
                data_in{1,ll}.regressor{1,25}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: complexity'};
                data_in{1,ll}.regressor{1,26}.name='compMAXinvFano12';     data_in{1,ll}.regressor{1,26}.value=[]; % SPE
                data_in{1,ll}.regressor{1,26}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: compMAXinvFano12'};
                data_in{1,ll}.regressor{1,27}.name='lowcompMAXinvFano12';     data_in{1,ll}.regressor{1,27}.value=[]; % SPE
                data_in{1,ll}.regressor{1,27}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: lowcompMAXinvFano12'};
                data_in{1,ll}.regressor{1,28}.name='highcompMAXinvFano12';     data_in{1,ll}.regressor{1,28}.value=[]; % SPE
                data_in{1,ll}.regressor{1,28}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: highcompMAXinvFano12'};
                data_in{1,ll}.regressor{1,29}.name='compinvFanoM1';     data_in{1,ll}.regressor{1,29}.value=[]; % SPE
                data_in{1,ll}.regressor{1,29}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: compinvFanoM1'};
                data_in{1,ll}.regressor{1,30}.name='compinvFanoM2';     data_in{1,ll}.regressor{1,30}.value=[]; % SPE
                data_in{1,ll}.regressor{1,30}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: compinvFanoM2'};
                data_in{1,ll}.regressor{1,31}.name='uncer';     data_in{1,ll}.regressor{1,31}.value=[]; % SPE
                data_in{1,ll}.regressor{1,31}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: uncer'};
                data_in{1,ll}.regressor{1,32}.name='medianMAXREL';     data_in{1,ll}.regressor{1,32}.value=[]; % SPE
                data_in{1,ll}.regressor{1,32}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: medianMAXREL'};
                data_in{1,ll}.regressor{1,33}.name='REL_COMP';     data_in{1,ll}.regressor{1,33}.value=[]; % SPE
                data_in{1,ll}.regressor{1,33}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: REL_COMP'};
                data_in{1,ll}.regressor{1,34}.name='REL_UNCER';     data_in{1,ll}.regressor{1,34}.value=[]; % SPE
                data_in{1,ll}.regressor{1,34}.Tag={'row1: session';'row2: block';'row3: trial in block';'row4: trial_s in trial';'row5: goal';'row6: active_model';'row7: REL_UNCER'};
                
                que_window=1;
                i_que=1;    NegLogLik_sarsa=zeros(1,que_window);    NegLogLik_fwd=zeros(1,que_window);
                data_in{1,ll}.model_error{1,1}.name={'Sum of negative logLikelihood of fwd model - the smaller the better fit'};  data_in{1,ll}.model_error{1,1}.value=[];
                data_in{1,ll}.model_error{1,2}.name={'Sum of negative logLikelihood of sarsa model - the smaller the better fit'};  data_in{1,ll}.model_error{1,2}.value=[];
                data_in{1,ll}.model_error{1,3}.name={'negative logLikelihood - the smaller the better fit'};  data_in{1,ll}.model_error{1,3}.value=[];
            end % #######################################
            
            
            
            
            Sum_NegLogLik_each_simul=0;
            Sum_NegLogLik_each_simul_eachMode=zeros(1,4);
            Num_occur_eachMode=zeros(1,4); % # of each mode (G,G',H,H')
            OBS=[];
            
            
            
            %% Simulation
            % (1) phase 1 - random action, no reward
            state_fwd=myState;  state_fwd.name='fwd';
            state_sarsa=myState;    state_sarsa.name='sarsa';
            state_sarsa_idi = myState; state_sarsa_idi.name='sarsa';
            
            if(data_in{1,ll}.map_type==1)
                map=myMap;  map0=myMap; map0_s=myMap;
            end
            if(data_in{1,ll}.map_type==2)
                map=myMap_new;  map0=myMap_new; map0_s=myMap_new;
            end
            
            
            %% (1) phase 1 - 'pre' training
            
            num_max_trial0=size(data_in{1,ll}.HIST_behavior_info_pre{1,1},1);
            map0.epoch=kk;                map0_s.epoch=kk;
            map0.data=data_in{1,ll}.HIST_behavior_info_pre{1,1};  map0_s.data=data_in{1,ll}.HIST_behavior_info_pre{1,1};
            
            opt_state_space.use_data=mode.experience_sbj_events(1); % use subject data for state-transition
            
            if(mode.experience_sbj_events(1)~=-1)
                
                state_sarsa.Qs=state_sarsa.Q;
                state_sarsa.Qr=state_sarsa.Q;
                state_sarsa.Qb=state_sarsa.Q;
                
                % fwd learning
                i=0;  cond=1;
                while ((i<num_max_trial0)&&(cond))
                    i=i+1;
                    map0.trial=i;
                    
                    Tprob_Ulow=[0.9 0.1];   Tprob_Uhigh=[0.5 0.5];
                    action_valid=[1 1 1 1]; % L1,R1,L2,R2
                    block_condition=data_in{1,ll}.HIST_block_condition_pre{1,1}(2,i);
                    rwd_condition=data_in{1,ll}.HIST_block_condition_pre{1,1}(3:5,i);
                    
                    switch block_condition
                        case 1
                            Tprob0=Tprob_Ulow;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            current_choice_type=2;
                        case 2
                            Tprob0=Tprob_Ulow;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            current_choice_type=1;
                        case 3
                            Tprob0=Tprob_Uhigh;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            current_choice_type=2;
                        case 4
                            Tprob0=Tprob_Uhigh;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            current_choice_type=1;
                    end
                    % update T_prob
                    map0=Model_Map_update(map0,'T',Tprob0);
                    % set Reward
                    map0=Model_Map_update(map0,'R',rwd_condition);
                    
                    
                    % initializing the state
                    param_fwd0=param_fwd;   param_fwd0.alpha=0.15;
                    [state_fwd map0]=StateClear(state_fwd,map0);
                    while (~state_fwd.JobComplete)
                        % decision
                        if(mode.experience_sbj_events(1)==1)
                            state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_behavior_data_save');
                        else
                            state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'decision_random');
                        end
                        % state transition
                        [state_fwd map0]=StateSpace_v1(state_fwd,map0,opt_state_space);  % map&state index ++
                        % 1. fwd model update
                        state_fwd=Model_RL(state_fwd, param_fwd0, map0, 'fwd_update');
                    end
                end
                
                % sarsa learning
                i=0;  cond=1;
                while ((i<num_max_trial0)&&(cond))
                    if mode.ori == 1
                        
                    else
                        
                        mode_data = map0.data(map0.trial,18:20);
                        [goal_color,goal_color_ind]  = max(mode_data);
                        state_sarsa.goal_color_ind = goal_color_ind;
                        silver = mode_data(1);
                        red = mode_data(2);
                        blue = mode_data(3);
                        
                        % Q switch for each colour
                        switch goal_color
                            case silver
                                state_sarsa.Q = state_sarsa.Qs;
                            case red
                                state_sarsa.Q = state_sarsa.Qr;
                            case blue
                                state_sarsa.Q = state_sarsa.Qb;
                        end
                        
                    end
                    
                    
                    i=i+1;
                    map0_s.trial=i;
                    
                    Tprob_Ulow=[0.9 0.1];   Tprob_Uhigh=[0.5 0.5];
                    action_valid=[1 1 1 1]; % L1,R1,L2,R2
                    block_condition=data_in{1,ll}.HIST_block_condition_pre{1,1}(2,i);
                    rwd_condition=data_in{1,ll}.HIST_block_condition_pre{1,1}(3:5,i);
                    
                    switch block_condition
                        case 1
                            Tprob0=Tprob_Ulow;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            current_choice_type=2;
                        case 2
                            Tprob0=Tprob_Ulow;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            current_choice_type=1;
                        case 3
                            Tprob0=Tprob_Uhigh;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            current_choice_type=2;
                        case 4
                            Tprob0=Tprob_Uhigh;
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            map0_s.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            current_choice_type=1;
                    end
                    % update T_prob
                    map0=Model_Map_update(map0,'T',Tprob0);
                    % set Reward
                    map0=Model_Map_update(map0,'R',rwd_condition);
                    
                    
                    
                    % initializing the state
                    param_sarsa0=param_sarsa;   param_sarsa0.alpha=param_fwd0.alpha*1.2; % make a bit faster than MB due to the slower convergence rate of MF.
                    [state_sarsa map0_s]=StateClear(state_sarsa,map0_s);
                    while (~state_sarsa.JobComplete)
                        % 0. current action selection : (s,a) - using arbitrator's Q-value
                        if(mode.experience_sbj_events(1)==1)
                            state_sarsa=Model_RL(state_sarsa, param_sarsa0, map0_s, 'decision_behavior_data_save');
                        else
                            state_sarsa=Model_RL(state_sarsa, param_sarsa0, map0_s, 'decision_random');
                        end
                        % 1. sarsa state update (get reward and next state) : (r,s')
                        [state_sarsa map0_s]=StateSpace_v1(state_sarsa,map0_s,opt_state_space); % map&state index ++
                        % 1. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                        state_sarsa=Model_RL(state_sarsa, param_sarsa0, map0_s, 'decision_hypo');
                        % 1. sarsa model upate
                        if mode.ori == 2
                            state_sarsa=Model_RL(state_sarsa, param_sarsa0, map0_s, 'sarsa_3q_update');
                        else
                            state_sarsa=Model_RL(state_sarsa, param_sarsa0, map0_s, 'sarsa_update');
                        end
                        
                        
                    end
                end
                
                % Q-value synchronization (subjects explored to learn about
                % the map. so use fwd Q as an initial Q-value
                %                 state_sarsa.Q=state_fwd.Q;
                
                % save current inital configuration
                data_in{1,ll}.init_state_fwd=state_fwd;
                data_in{1,ll}.init_state_sarsa=state_sarsa;
                
            else % =-1: retrieve saved configuration
                state_fwd.Q=data_in{1,ll}.init_state_fwd.Q;
                state_sarsa.Q=data_in{1,ll}.init_state_sarsa.Q;
            end
            
            
            
            %                 disp('- pretraining completed.')
            
            %% (2) phase 2 - intact action, rewards given
            
            %             param_sarsa.alpha=param_sarsa.alpha/0.5;
            %             param_fwd.alpha=param_fwd.alpha/0.5;
            
            num_max_session=size(data_in{1,ll}.HIST_behavior_info,2);
            
            
            cond=1;
            myArbitrator_top=myArbitrator;
            mode_data_prev=[5 5 5];
            block_cond_prev=-1; on_h=-1; on_g=-1;
            
            % initialize Q for each color
            state_sarsa.Qs=state_sarsa.Q;
            state_sarsa.Qr=state_sarsa.Q;
            state_sarsa.Qb=state_sarsa.Q;
            
            
            for ind_sess=1:1:num_max_session
                % enter each session data into map
                i=0;
                num_max_trial=size(data_in{1,ll}.HIST_behavior_info{1,ind_sess},1);
                map.epoch=kk;
                map.data=data_in{1,ll}.HIST_behavior_info{1,ind_sess};
                map.condition_id=data_in{1,ll}.HIST_block_condition{1,ind_sess}(2,:);
                
                
                
                %                     disp(sprintf('- Sbj[%d], Simul[%d/%d], Session [%d/%d]...',###,kk,mode.total_simul,ind_sess,num_max_session));
                
                while ((i<num_max_trial)&&(cond))
                    i=i+1;
                    map.trial=i;
                    %                         disp(sprintf('- Simul[%d/%d], Session [%d/%d], Trial [%d/%d]...',kk,mode.total_simul,ind_sess,num_max_session,i,num_max_trial));
                    
                    % initializing the state
                    [state_fwd map]=StateClear(state_fwd,map);
                    [state_sarsa map]=StateClear(state_sarsa,map);
                    
                    % read a reward value set from data : 1x3
                    mode_data=map.data(map.trial,18:20);
                    
                    % color define
                    silver = map.data(map.trial,18);
                    red = map.data(map.trial,19);
                    blue = map.data(map.trial,20);
                    [goal_color,goal_color_ind]  = max(mode_data);
                    state_sarsa.goal_color_ind = goal_color_ind;
                    
                    state_2 = map.data(map.trial,5);
                    gc_ct = goal_color;
                    if i ==1
                        gc_pt = goal_color + 1 ;% just make it different from the gc_ct
                    end
                    % Q switch for each colour
                    if mode.ori == 1
                        
                    else
                        switch goal_color
                            case silver
                                state_sarsa.Q = state_sarsa.Qs;
                            case red
                                state_sarsa.Q = state_sarsa.Qr;
                            case blue
                                state_sarsa.Q = state_sarsa.Qb;
                        end
                    end
                    
                    % read block and rwd value set
                    block_condition=data_in{1,ll}.HIST_block_condition{1,ind_sess}(2,map.trial);
                    
                    % rvsd complexity or not
                    if mode.rvs == 1
                        high_comp = 1;
                        low_comp = 2;
                    else
                        high_comp = 2;
                        low_comp = 1;
                    end
                    
                    
                    
                    % get available action set
                    switch block_condition
                        case 1
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            myArbitrator_top.var_cogload=low_comp;
                        case 2
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            myArbitrator_top.var_cogload=high_comp;
                        case 3
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 -Inf -Inf]];
                            myArbitrator_top.var_cogload=low_comp;
                        case 4
                            map.action_valid=[[1 1 -Inf -Inf]; [1 1 1 1]];
                            myArbitrator_top.var_cogload=high_comp;
                    end
                    
                    % mode change whenever there is a visual cue change
                    case_number=0;
                    myArbitrator_top.backward_flag=0;
                    if(sum(mode_data_prev~=mode_data)>0) % if there is any visual cue change, then goal mode
                        myArbitrator_top.backward_flag=1; % after backward update, should set it to 0.
                        %                         myArbitrator_top.ind_active_model=1; % switching the mode
                        %                         myArbitrator_top.m1_prob_prev=0.9; % changing the choice prob accordingly
                        %                         myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        %                         if((mode_data_prev~=-1)&&(mode_data==-1)) % goal mode -> habitual mode
                        %                             myArbitrator_top.ind_active_model=2; % switching the mode
                        %                             myArbitrator_top.m1_prob_prev=0.1; % changing the choice prob accordingly
                        %                             myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        %                         end
                    end
                    if(mode.USE_FWDSARSA_ONLY==1) % forward only
                        myArbitrator_top.ind_active_model=1; % switching the mode
                        myArbitrator_top.m1_prob_prev=0.99; % changing the choice prob accordingly
                        myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        myArbitrator_top.Time_Step=1e-20; % extremely slow, so do not switch to the other learner
                    end
                    if(mode.USE_FWDSARSA_ONLY==2) % sarsa only
                        myArbitrator_top.ind_active_model=2; % switching the mode
                        myArbitrator_top.m1_prob_prev=0.01; % changing the choice prob accordingly
                        myArbitrator_top.m1_prob=myArbitrator_top.m1_prob_prev;
                        myArbitrator_top.Time_Step=1e-20; % extremely slow, so do not switch to the other learner
                    end
                    myArbitrator_top.m2_prob=1-myArbitrator_top.m1_prob;
                    
                    if(mode.out==0) % debug
                        OBS=[OBS [myArbitrator_top.ind_active_model; mode_data'; myArbitrator_top.m1_prob; state_fwd.Q(1,1); state_fwd.Q(1,2); state_sarsa.Q(1,1); state_sarsa.Q(1,2)]];
                    end
                    %
                    
                    opt_state_space.use_data=mode.experience_sbj_events(2); % use subject data for state-transition
                    
                    while (((myArbitrator_top.ind_active_model==1)&&(~map.JobComplete))||((myArbitrator_top.ind_active_model==2)&&(~map.JobComplete)))
                        
                        
                        % index synchronization
                        state_fwd.index=map.index;                state_sarsa.index=map.index;
                        
                        %% fwd mode: backward update of fwd model
                        if(myArbitrator_top.backward_flag==1) % if the agent detects context change
                            % (1) revaluation : reevaluation for goal mode
                            map=Model_Map_update(map,'R',mode_data');
                            % (2) backward update of the fwd model
                            if(mode.USE_BWDupdate_of_FWDmodel==1)
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'bwd_update');
                            end
                        else
                            %
                        end
                        
                        
                        % regressor part ############################
                        % DO NOT change the order !!!!!!!!!!!!!!!
                        % trial_s in trial (1 only)
                        % regressor "right before" making the first choice at the first state
                        
                        if(mode.out==99) % regressors
                            if(map.index==1) % only for the *first state*
                                an_s=ind_sess; % session
                                an_i=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,1);% block in session
                                an_ib=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,2);% trial in block
                                an_g=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,18); % goal state
                                % regressor1
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; SPE
                                data_in{1,ll}.regressor{1,1}.value=[data_in{1,ll}.regressor{1,1}.value vec_reg];
                                % regressor2
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; RPE
                                data_in{1,ll}.regressor{1,2}.value=[data_in{1,ll}.regressor{1,2}.value vec_reg];
                                % regressor3
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_var;]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m1 (3x1);
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m1 (3x1);
                                data_in{1,ll}.regressor{1,3}.value=[data_in{1,ll}.regressor{1,3}.value vec_reg];
                                % regressor4
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_var;]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m2 (3x1);
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m2 (3x1);
                                data_in{1,ll}.regressor{1,4}.value=[data_in{1,ll}.regressor{1,4}.value vec_reg];
                                % regressor5
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_mean;]; % session; block; trial; trial_s; goal; active_model; mean of m1 (3x1);
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; mean of m1 (3x1);
                                data_in{1,ll}.regressor{1,5}.value=[data_in{1,ll}.regressor{1,5}.value vec_reg];
                                % regressor6
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_mean;]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.regressor{1,6}.value=[data_in{1,ll}.regressor{1,6}.value vec_reg];
                                % regressor7
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_inv_Fano;]; % session; block; trial; trial_s; goal; active_model; invFano of m1 (3x1);
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; invFano of m1 (3x1);
                                data_in{1,ll}.regressor{1,7}.value=[data_in{1,ll}.regressor{1,7}.value vec_reg];
                                % regressor8
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; invFano of m2 (3x1);
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_inv_Fano;]; % session; block; trial; trial_s; goal; active_model; invFano of m2 (3x1);
                                data_in{1,ll}.regressor{1,8}.value=[data_in{1,ll}.regressor{1,8}.value vec_reg];
                                % regressor9
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_prob;]; % session; block; trial; trial_s; goal; active_model; choice prob of m1
                                data_in{1,ll}.regressor{1,9}.value=[data_in{1,ll}.regressor{1,9}.value vec_reg];
                                % regressor10
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_prob;]; % session; block; trial; trial_s; goal; active_model; choice prob of m2
                                data_in{1,ll}.regressor{1,10}.value=[data_in{1,ll}.regressor{1,10}.value vec_reg];
                                % regressor11 - Qfwd
                                an_j=1;
                                action_chosen=map.data(map.trial,6+1);
                                if(action_chosen==1)    action_unchosen=2;  end
                                if(action_chosen==2)    action_unchosen=1;  end
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_fwd.Q(1,action_chosen)-state_fwd.Q(1,action_unchosen)]; % session; block; trial; trial_s; goal; active_model; x
                                %                                 [tmp_max action_chosen_fwd]=max(state_fwd.Q(1,:));  [tmp_min action_unchosen_fwd]=min(state_fwd.Q(1,:));
                                action_chosen_fwd=action_chosen;    mul_val=1;%myArbitrator_top.m1_prob;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; mul_val*state_fwd.Q(1,action_chosen_fwd)]; % session; block; trial; trial_s; goal; active_model; x
                                data_in{1,ll}.regressor{1,11}.value=[data_in{1,ll}.regressor{1,11}.value vec_reg];
                                % regressor12 - Qsarsa
                                an_j=1;
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_sarsa.Q(1,action_chosen)-state_sarsa.Q(1,action_unchosen)]; % session; block; trial; trial_s; goal; active_model; x
                                %                                 [tmp_max action_chosen_sarsa]=max(state_sarsa.Q(1,:));  [tmp_min action_unchosen_sarsa]=min(state_sarsa.Q(1,:));
                                action_chosen_sarsa=action_chosen;      mul_val=1;%myArbitrator_top.m2_prob;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; mul_val*state_sarsa.Q(1,action_chosen_sarsa)]; % session; block; trial; trial_s; goal; active_model; x
                                data_in{1,ll}.regressor{1,12}.value=[data_in{1,ll}.regressor{1,12}.value vec_reg];
                                % regressor13 - Qarb
                                an_j=1;
                                myArbitrator_top.Q=...
                                    ((myArbitrator_top.m1_prob*state_fwd.Q).^myArbitrator_top.p+...
                                    (myArbitrator_top.m2_prob*state_sarsa.Q).^myArbitrator_top.p).^(1/myArbitrator_top.p);
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.Q(1,action_chosen)-myArbitrator_top.Q(1,action_unchosen)]; % session; block; trial; trial_s; goal; active_model; choice prob of m1
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.Q(1,action_chosen)]; % session; block; trial; trial_s; goal; active_model; choice prob of m1
                                data_in{1,ll}.regressor{1,13}.value=[data_in{1,ll}.regressor{1,13}.value vec_reg];
                                % regressor14 - dQbwdEnergy (norm of Q-vector change)
                                an_j=1;
                                if(myArbitrator_top.backward_flag==1) % if the agent detects context change > bwd_updated
                                    vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_fwd.dQ_bwd_energy]; % session; block; trial; trial_s; goal; active_model; x
                                else
                                    vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                                end
                                data_in{1,ll}.regressor{1,14}.value=[data_in{1,ll}.regressor{1,14}.value vec_reg];
                                % regressor15 - dQbwd (mean)
                                an_j=1;
                                if(myArbitrator_top.backward_flag==1) % if the agent detects context change > bwd_updated
                                    vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_fwd.dQ_mean_energy]; % session; block; trial; trial_s; goal; active_model; x
                                else
                                    vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                                end
                                data_in{1,ll}.regressor{1,15}.value=[data_in{1,ll}.regressor{1,15}.value vec_reg];
                                % regressor16 - duncertaintyM1 (change of M1 uncertainty)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; x
                                data_in{1,ll}.regressor{1,16}.value=[data_in{1,ll}.regressor{1,16}.value vec_reg];
                                % regressor17 - duncertaintyM1 (change of M1 invFano)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; x
                                data_in{1,ll}.regressor{1,17}.value=[data_in{1,ll}.regressor{1,17}.value vec_reg];
                                % regressor18 - transition rate for MF2MB(alpha)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.transition_rate21]; % session; block; trial; trial_s; goal; active_model; alpha (MF->MB) based on the performance of MF
                                data_in{1,ll}.regressor{1,18}.value=[data_in{1,ll}.regressor{1,18}.value vec_reg];
                                % regressor19 - transition rate for MB2MF(beta)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.transition_rate12]; % session; block; trial; trial_s; goal; active_model; % beta (MB->MF) based on the performance of MB
                                data_in{1,ll}.regressor{1,19}.value=[data_in{1,ll}.regressor{1,19}.value vec_reg];
                                % regressor20- (invF1-invF2)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; invFano of m1 - m2 (1x1);
                                data_in{1,ll}.regressor{1,20}.value=[data_in{1,ll}.regressor{1,20}.value vec_reg];
                                % regressor21- abs(invF1-invF2)
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; zeros(3,1)]; % session; block; trial; trial_s; goal; active_model; abs(invFano of m1 - m2) (1x1);
                                data_in{1,ll}.regressor{1,21}.value=[data_in{1,ll}.regressor{1,21}.value vec_reg];
                                % regressor22- max(invF1-invF2)
                                % will be added at the end of this code
                                % regressor23- conflict btw invF1&invF2
                                % will be added at the end of this code
                                
                                % regressor24
                                an_j=1;
                                %                                 len1 = sqrt(mode_data_prev(1)^2 + mode_data_prev(2)^2 + mode_data_prev(3)^2);
                                %                                 len2 = sqrt(mode_data(1)^2 + mode_data(2)^2 + mode_data(3)^2);
                                %                                 result = 1 / (sum(mode_data_prev.*mode_data) ./ (len1 * len2));
                                %                                 result = result - 1;
                                %                                 gc_consist = result;
                                if gc_pt == gc_ct
                                    gc_consist = 1;
                                else
                                    gc_consist = 0;
                                    
                                end
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; gc_consist]; % session; block; trial; trial_s; goal; active_model; SPE
                                data_in{1,ll}.regressor{1,24}.value=[data_in{1,ll}.regressor{1,24}.value vec_reg];
                                
                                
                                % regressor25- complexity
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; abs(invFano of m1 - m2) (1x1);
                                data_in{1,ll}.regressor{1,25}.value=[data_in{1,ll}.regressor{1,25}.value vec_reg];
                            end
                            
                        end
                        
                        
                        
                        %% Compute negative log-likelihood : evaluate the arbitrator softmax using sbj's state,action
                        state_data=map.data(map.trial,3+map.index);
                        action_chosen=map.data(map.trial,6+map.index);
                        % compute real Q-value by merging two Q (fwd & sarsa)
                        myArbitrator_top.Q=...
                            ((myArbitrator_top.m1_prob*state_fwd.Q).^myArbitrator_top.p+...
                            (myArbitrator_top.m2_prob*state_sarsa.Q).^myArbitrator_top.p).^(1/myArbitrator_top.p);
                        myArbitrator_top.Q_old=myArbitrator_top.Q;
                        q_in=(map.action_valid(map.index,:)).*myArbitrator_top.Q(state_data,:);
                        % exploitative
                        %                         var_exp=exp(myArbitrator_top.var_cogload*myArbitrator_top.tau_softmax*q_in); % (N_actionx1)
                        % explorative
                        var_exp=exp((myArbitrator_top.var_cogload.^(-1) *2)*myArbitrator_top.tau_softmax*q_in); % (N_actionx1)
                        eval_num=log(var_exp(action_chosen)/sum(var_exp));
                        Sum_NegLogLik=Sum_NegLogLik-eval_num/(mode.total_simul);
                        Sum_NegLogLik_each_simul=Sum_NegLogLik_each_simul-eval_num/(1*1);
                        % block condition of each trial: 1:G(low T uncertainty), 2:G'(high T uncertainty), 3:H(high T uncertainty), 4:H'(low T uncertainty)
                        block_cond=data_in{1,ll}.HIST_block_condition{1,ind_sess}(2,map.trial);
                        Sum_NegLogLik_each_simul_eachMode(block_cond)=Sum_NegLogLik_each_simul_eachMode(block_cond)-eval_num/(mode.total_simul*tot_num_sbj);
                        Num_occur_eachMode(block_cond)=Num_occur_eachMode(block_cond)+1;
                        
                        
                        %% compute negative log-likelihood in time window for fwd&sarsa separately
                        if(mode.out==99) % in regressor generation mode, do simulation only once
                            if(map.index==1) % only for the *first state*
                                i_que=i_que+1;                            i_que=min(que_window,i_que+1);
                                %                             if(i_que<que_window)
                                %                             else
                                state_data=map.data(map.trial,3+map.index);
                                action_chosen=map.data(map.trial,6+1);
                                % compute real Q-value by merging two Q (fwd & sarsa)
                                q_in=(map.action_valid(map.index,:)).*state_fwd.Q(state_data,:);
                                var_exp_fwd=exp(param_fwd.tau*q_in); % (N_actionx1)
                                eval_num_fwd=log(var_exp_fwd(action_chosen)/sum(var_exp_fwd));
                                if(que_window>1)                                    NegLogLik_fwd(1:end-1)=NegLogLik_fwd(2:end);                                end
                                NegLogLik_fwd(end)=(-1.0)*eval_num_fwd;
                                q_in=(map.action_valid(map.index,:)).*state_sarsa.Q(state_data,:);
                                var_exp_sarsa=exp(param_sarsa.tau*q_in); % (N_actionx1)
                                eval_num_sarsa=log(var_exp_sarsa(action_chosen)/sum(var_exp_sarsa));
                                if(que_window>1)                                    NegLogLik_sarsa(1:end-1)=NegLogLik_sarsa(2:end);                                end
                                NegLogLik_sarsa(end)=(-1.0)*eval_num_sarsa;
                                %
                                q_in=(map.action_valid(map.index,:)).*myArbitrator_top.Q(state_data,:);
                                var_exp=exp(myArbitrator_top.tau_softmax*q_in); % (N_actionx1)
                                eval_num=log(var_exp(action_chosen)/sum(var_exp));
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; sum(NegLogLik_fwd)]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,1}.value=[data_in{1,ll}.model_error{1,1}.value vec_reg];
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; sum(NegLogLik_sarsa)]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,2}.value=[data_in{1,ll}.model_error{1,2}.value vec_reg];
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; eval_num]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,3}.value=[data_in{1,ll}.model_error{1,3}.value vec_reg];
                                %                             end
                            end
                        end
                        
                        
                        %% main computation
                        
                        % controlled update = simultaneous update of both model
                        QQ_prev=state_fwd.Q;
                        if(myArbitrator_top.ind_active_model==1) % fwd
                            % 0. current action selection : (s,a) - using arbitrator's Q-value
                            if(mode.experience_sbj_events(2)==1)
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_behavior_data_save');
                            else
                                state_fwd=Model_RL(state_fwd, param_fwd, map, 'decision_arbitrator', myArbitrator_top, state_sarsa);
                            end
                            
                            % 1. fwd state update (get reward and next state) : (r,s')
                            [state_fwd map]=StateSpace_v1(state_fwd,map,opt_state_space); % map&state index ++
                            state_sarsa.index=state_fwd.index;
                            % 1. fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
                            % 2. state synchronization
                            state_sarsa.state_history(state_fwd.index)=state_fwd.state_history(state_fwd.index);
                            state_sarsa.SARSA=state_fwd.SARSA;
                            % 3. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                            % 3. sarsa model upate
                            if mode.ori == 2
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_3q_update');
                            elseif mode.ori ==3
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_3q_update');
                                state_sarsa_idi=Model_RL(state_sarsa_idi, param_sarsa, map, 'sarsa_update');
                            else
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
                            end
                            % history synchronization
                            myArbitrator_top.state_history=state_fwd.state_history;
                            myArbitrator_top.reward_history=state_fwd.reward_history;
                        end
                        QQ_current=state_fwd.Q;
                        
                        if(myArbitrator_top.ind_active_model==2) % sarsa
                            % 0. current action selection : (s,a) - using arbitrator's Q-value
                            if(mode.experience_sbj_events(2)==1)
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_behavior_data_save');
                            else
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_arbitrator', myArbitrator_top, state_fwd);
                            end
                            % 1. sarsa state update (get reward and next state) : (r,s')
                            [state_sarsa map]=StateSpace_v1(state_sarsa,map,opt_state_space); % map&state index ++
                            state_fwd.index=state_sarsa.index;
                            % 1. sarsa next action selection : (s',a') - if s' is terminal, then no decision
                            state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'decision_hypo');
                            % 1. sarsa model upate
                            if mode.ori == 2
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_3q_update');
                            else
                                state_sarsa=Model_RL(state_sarsa, param_sarsa, map, 'sarsa_update');
                            end                            % 2. state synchronization
                            state_fwd.state_history(state_sarsa.index)=state_sarsa.state_history(state_sarsa.index);
                            state_fwd.SARSA=state_sarsa.SARSA;
                            % 3. fwd model update
                            state_fwd=Model_RL(state_fwd, param_fwd, map, 'fwd_update');
                            % history synchronization
                            myArbitrator_top.state_history=state_sarsa.state_history;
                            myArbitrator_top.reward_history=state_sarsa.reward_history;
                        end
                        
                        % ARBITRATOR: transition rate of the sarsa only
                        %                         [myArbitrator_top, state_fwd, state_sarsa]=Bayesian_Arbitration_v3(myArbitrator_top, state_fwd, state_sarsa, map); % full bayesian
                        [myArbitrator_top, state_fwd, state_sarsa]=Bayesian_Arbitration_v7(myArbitrator_top, state_fwd, state_sarsa, map, mode); % bayesian for m1, RPE estimator for m2 + cogload
                        
                        
                        % regressor part ############################
                        % DO NOT change the order !!!!!!!!!!!!!!!
                        % trial_s in trial (2 or 3)
                        % regressor "right AFTER" making a choice - so the "active model" is the model choice after the event
                        if(mode.out==99) % regressors
                            an_s=ind_sess; % session
                            an_i=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,1);% block in session
                            an_ib=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,2);% trial in block
                            an_g=data_in{1,ll}.HIST_behavior_info{1,ind_sess}(map.trial,18); % goal state
                            % regressor1
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_fwd.SPE_history(map.index)]; % session; block; trial; trial_s; goal; active_model; SPE
                            data_in{1,ll}.regressor{1,1}.value=[data_in{1,ll}.regressor{1,1}.value vec_reg];
                            % regressor2
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            if mode.ori == 2
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_sarsa.RPE_history(map.index, state_sarsa.goal_color_ind)]; % session; block; trial; trial_s; goal; active_model; RPE
                            else
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_sarsa.RPE_history(map.index)]; % session; block; trial; trial_s; goal; active_model; RPE
                            end
                            data_in{1,ll}.regressor{1,2}.value=[data_in{1,ll}.regressor{1,2}.value vec_reg];
                            % regressor3
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_var;]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m1 (3x1);
                            data_in{1,ll}.regressor{1,3}.value=[data_in{1,ll}.regressor{1,3}.value vec_reg];
                            % regressor4
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_var;]; % session; block; trial; trial_s; goal; active_model; Uncertaint of m2 (3x1);
                            data_in{1,ll}.regressor{1,4}.value=[data_in{1,ll}.regressor{1,4}.value vec_reg];
                            % regressor5
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_mean;]; % session; block; trial; trial_s; goal; active_model; mean of m1 (3x1);
                            data_in{1,ll}.regressor{1,5}.value=[data_in{1,ll}.regressor{1,5}.value vec_reg];
                            % regressor6
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_mean;]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                            data_in{1,ll}.regressor{1,6}.value=[data_in{1,ll}.regressor{1,6}.value vec_reg];
                            % regressor7
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_inv_Fano;]; % session; block; trial; trial_s; goal; active_model; invFano of m1 (3x1);
                            data_in{1,ll}.regressor{1,7}.value=[data_in{1,ll}.regressor{1,7}.value vec_reg];
                            % regressor8
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_inv_Fano;]; % session; block; trial; trial_s; goal; active_model; invFano of m2 (3x1);
                            data_in{1,ll}.regressor{1,8}.value=[data_in{1,ll}.regressor{1,8}.value vec_reg];
                            % regressor9
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m1_prob;]; % session; block; trial; trial_s; goal; active_model; choice prob of m1
                            data_in{1,ll}.regressor{1,9}.value=[data_in{1,ll}.regressor{1,9}.value vec_reg];
                            % regressor10
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.m2_prob;]; % session; block; trial; trial_s; goal; active_model; choice prob of m2
                            data_in{1,ll}.regressor{1,10}.value=[data_in{1,ll}.regressor{1,10}.value vec_reg];
                            % regressor11 - Qfwd
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            if(an_j==2)
                                state_data=map.data(map.trial,3+map.index);                                action_chosen=map.data(map.trial,6+2);
                                if(action_chosen==1)    action_unchosen=2;  end
                                if(action_chosen==2)    action_unchosen=1;  end
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_fwd.Q(state_data,action_chosen)-state_fwd.Q(state_data,action_unchosen)]; % session; block; trial; trial_s; goal; active_model; x
                                %                                 [tmp_max action_chosen_fwd]=max(state_fwd.Q(state_data,:));  [tmp_min action_unchosen_fwd]=min(state_fwd.Q(state_data,:));
                                action_chosen_fwd=action_chosen;        mul_val=1;%myArbitrator_top.m1_prob;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; mul_val*state_fwd.Q(state_data,action_chosen_fwd)]; % session; block; trial; trial_s; goal; active_model; x
                            else % an_j==3: no Q-value because it is the last state
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                            end
                            data_in{1,ll}.regressor{1,11}.value=[data_in{1,ll}.regressor{1,11}.value vec_reg];
                            % regressor12 - Qsarsa
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            if(an_j==2)
                                state_data=map.data(map.trial,3+map.index);                                action_chosen=map.data(map.trial,6+2);
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; state_sarsa.Q(state_data,action_chosen)-state_sarsa.Q(state_data,action_unchosen)]; % session; block; trial; trial_s; goal; active_model; x
                                %                                 [tmp_max action_chosen_sarsa]=max(state_sarsa.Q(state_data,:));  [tmp_min action_unchosen_sarsa]=min(state_sarsa.Q(state_data,:));
                                action_chosen_sarsa=action_chosen;      mul_val=1;%myArbitrator_top.m2_prob;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; mul_val*state_sarsa.Q(state_data,action_chosen_sarsa)]; % session; block; trial; trial_s; goal; active_model; x
                            else % an_j==3: no Q-value because it is the last state
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                            end
                            data_in{1,ll}.regressor{1,12}.value=[data_in{1,ll}.regressor{1,12}.value vec_reg];
                            % regressor13 - Qarb
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            if(an_j==2)
                                if map.action_valid(2,3) == -inf
                                    action_unchosen_arb = action_unchosen;
                                else
                                    switch action_chosen
                                        case 1
                                            action_unchosen_arb = 2;
                                        case 2
                                            action_unchosen_arb = 1;
                                        case 3
                                            action_unchosen_arb = 4;
                                        case 4
                                            action_unchosen_arb = 3;
                                    end
                                end
                                %
                                %                                 if map.action_valid(2,3) == -inf
                                %                                     action_unchosen_arb = action_unchosen;
                                %                                 else
                                %                                     action_unchosen_arb = setdiff([1,2,3,4],action_chosen);
                                %
                                %                                 end
                                state_data=map.data(map.trial,3+map.index);                                action_chosen=map.data(map.trial,6+2);
                                %                                 myArbitrator_top.Q=...
                                %                                     ((myArbitrator_top.m1_prob*state_fwd.Q).^myArbitrator_top.p+...
                                %                                     (myArbitrator_top.m2_prob*state_sarsa.Q).^myArbitrator_top.p).^(1/myArbitrator_top.p);
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.Q(state_data,action_chosen)]; % session; block; trial; trial_s; goal; active_model; x
                                %                                 vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.Q(state_data,action_chosen)-myArbitrator_top.Q(state_data,action_unchosen_arb)]; % session; block; trial; trial_s; goal; active_model; x
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.Q(state_data,action_chosen)-myArbitrator_top.Q(state_data,action_unchosen_arb)]; % session; block; trial; trial_s; goal; active_model; x
                                
                            else % an_j==3: no Q-value because it is the last state
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                            end
                            data_in{1,ll}.regressor{1,13}.value=[data_in{1,ll}.regressor{1,13}.value vec_reg];
                            % regressor14 - dQbwdEnergy (norm of Q-vector change) -> no bwd update
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                            data_in{1,ll}.regressor{1,14}.value=[data_in{1,ll}.regressor{1,14}.value vec_reg];
                            % regressor15 - dQbwd (mean) -> no bwd update
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; x
                            data_in{1,ll}.regressor{1,15}.value=[data_in{1,ll}.regressor{1,15}.value vec_reg];
                            % regressor16 - duncertaintyM1 (change of M1 uncertainty)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; (myArbitrator_top.m1_var-myArbitrator_top.m1_var_old)];
                            data_in{1,ll}.regressor{1,16}.value=[data_in{1,ll}.regressor{1,16}.value vec_reg];
                            % regressor17 - duncertaintyM1 (change of M1 invFano)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; (myArbitrator_top.m1_inv_Fano-myArbitrator_top.m1_inv_Fano_old)];
                            data_in{1,ll}.regressor{1,17}.value=[data_in{1,ll}.regressor{1,17}.value vec_reg];
                            % regressor18 - transition rate for MF2MB(alpha)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.transition_rate21]; % session; block; trial; trial_s; goal; active_model; alpha (MF->MB) based on the performance of MF
                            data_in{1,ll}.regressor{1,18}.value=[data_in{1,ll}.regressor{1,18}.value vec_reg];
                            % regressor19 - transition rate for MB2MF(beta)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; myArbitrator_top.transition_rate12]; % session; block; trial; trial_s; goal; active_model; % beta (MB->MF) based on the performance of MB
                            data_in{1,ll}.regressor{1,19}.value=[data_in{1,ll}.regressor{1,19}.value vec_reg];
                            % regressor20- (invF1-invF2)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            input0=myArbitrator_top.m1_inv_Fano;    input1=myArbitrator_top.m1_wgt'*myArbitrator_top.m1_inv_Fano;   invF1normalized=input1/sum(input0);
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; invF1normalized-myArbitrator_top.m2_inv_Fano]; % session; block; trial; trial_s; goal; active_model; invFano of m1 - m2 (3x1); [note]invF2 is already normalized
                            data_in{1,ll}.regressor{1,20}.value=[data_in{1,ll}.regressor{1,20}.value vec_reg];
                            % regressor21- abs(invF1-invF2)
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; abs(invF1normalized-myArbitrator_top.m2_inv_Fano)]; % session; block; trial; trial_s; goal; active_model; abs of invFano of m1 - m2 (3x1);
                            data_in{1,ll}.regressor{1,21}.value=[data_in{1,ll}.regressor{1,21}.value vec_reg];
                            % regressor22- max(invF1-invF2)
                            % will be added at the end of this code
                            % regressor23- conflict btw invF1&invF2
                            % will be added at the end of this code
                            
                            % regressor24
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; 0]; % session; block; trial; trial_s; goal; active_model; SPE
                            data_in{1,ll}.regressor{1,24}.value=[data_in{1,ll}.regressor{1,24}.value vec_reg];
                            
                            % regressor25
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; (myArbitrator_top.var_cogload)*1]; % session; block; trial; trial_s; goal; active_model; SPE
                            data_in{1,ll}.regressor{1,25}.value=[data_in{1,ll}.regressor{1,25}.value vec_reg];
                            
                            
                        end
                        
                        %% compute negative log-likelihood in time window for fwd&sarsa separately
                        if(mode.out==99) % in regressor generation mode, do simulation only once
                            state_data=map.data(map.trial,3+map.index);
                            an_j=myArbitrator_top.index+1; % trial_s in trial (2 or 3)
                            action_chosen=map.data(map.trial,6+2);
                            if(an_j==3) % do not compute LogLik at the last block : copy from the previous block - assuming that subject will be remain with the same model
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; NaN]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,1}.value=[data_in{1,ll}.model_error{1,1}.value vec_reg];
                                data_in{1,ll}.model_error{1,2}.value=[data_in{1,ll}.model_error{1,2}.value vec_reg];
                                data_in{1,ll}.model_error{1,3}.value=[data_in{1,ll}.model_error{1,3}.value vec_reg];
                            else % an_j==2
                                i_que=i_que+1;                            i_que=min(que_window,i_que+1);
                                % compute real Q-value by merging two Q (fwd & sarsa)
                                q_in=(map.action_valid(map.index,:)).*state_fwd.Q(state_data,:);
                                var_exp_fwd=exp(param_fwd.tau*q_in); % (N_actionx1)
                                eval_num_fwd=log(var_exp_fwd(action_chosen)/sum(var_exp_fwd));
                                if(que_window>1)                                    NegLogLik_fwd(1:end-1)=NegLogLik_fwd(2:end);                                end
                                NegLogLik_fwd(end)=(-1.0)*eval_num_fwd;
                                q_in=(map.action_valid(map.index,:)).*state_sarsa.Q(state_data,:);
                                var_exp_sarsa=exp(param_sarsa.tau*q_in); % (N_actionx1)
                                eval_num_sarsa=log(var_exp_sarsa(action_chosen)/sum(var_exp_sarsa));
                                if(que_window>1)                                NegLogLik_sarsa(1:end-1)=NegLogLik_sarsa(2:end);                                end
                                NegLogLik_sarsa(end)=(-1.0)*eval_num_sarsa;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; sum(NegLogLik_fwd)]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,1}.value=[data_in{1,ll}.model_error{1,1}.value vec_reg];
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; sum(NegLogLik_sarsa)]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,2}.value=[data_in{1,ll}.model_error{1,2}.value vec_reg];
                                q_in=(map.action_valid(map.index,:)).*myArbitrator_top.Q(state_data,:);
                                var_exp=exp(myArbitrator_top.tau_softmax*q_in); % (N_actionx1)
                                eval_num=log(var_exp(action_chosen)/sum(var_exp));
                                an_j=1;
                                vec_reg=[an_s; an_i; an_ib; an_j; an_g; myArbitrator_top.ind_active_model; eval_num]; % session; block; trial; trial_s; goal; active_model; mean of m2 (3x1);
                                data_in{1,ll}.model_error{1,3}.value=[data_in{1,ll}.model_error{1,3}.value vec_reg];
                            end
                        end
                        
                        % DISPLAY Q-VAL change according to goal change
                        if(mode.DEBUG_Q_VALUE_CHG==1)
                            disp(sprintf('goal: %d->%d, Q_debug=[state|| fwdQ(L)|fwdQ(R) -> fwdQ(L)|fwdQ(R) | arbQ(L)|arbQ(R)]',mode_data_prev, mode_data))
                            %                                 Q_debug=[[1:1:9]' QQ_prev QQ_current myArbitrator_top.Q(:,:)]
                            %                                 myArbitrator_top.m1_prob
                            %                                 state_fwd.SARSA
                            [state_data, action_chosen, var_exp, (-1)*eval_num]
                            input('pause...')
                        end
                        
                        
                    end
                    
                    % save to previous mode data
                    mode_data_prev=mode_data;
                    
                    % update Q value for each color
                    %                     if (myArbitrator_top.var_cogload==1) && ((goal_color==blue) || (goal_color==red && state_2==2))
                    %                             dummy_q = state_sarsa.Q;
                    %                     else
                    switch goal_color
                        case silver
                            state_sarsa.Qs = state_sarsa.Q;
                        case red
                            state_sarsa.Qr = state_sarsa.Q;
                        case blue
                            state_sarsa.Qb = state_sarsa.Q;
                    end
                    %                     end
                    gc_pt = goal_color;
                end
                
                
            end
            
            if(mode.simul_process_display==1)
                % display
                %             msg=[sprintf('- processing subject(%d/%d)',ll,tot_num_sbj) repmat('.',1,kk)];
                msg=[sprintf('       - processing subject(%d/%d),simulation(%d/%d),NegLogLik=%04.2f',ll,tot_num_sbj,kk,mode.total_simul,Sum_NegLogLik_each_simul)];
                fprintf(repmat('\b',1,n_tmp));            fprintf(msg);            n_tmp=numel(msg);
            end
            
            
        end
        if(mode.simul_process_display==1)        fprintf('\n');     end
    end
    
    %% range/mean-correction for diff invF regressors
    if(mode.out==99) % regressors
        for ll=1:1:tot_num_sbj % each subject
            zero_pt_ind=find(mod([1:1:length(data_in{1,ll}.regressor{1,7}.value(8,:))]+2,3)==0);
            % get max, min
            rr_min1=min(data_in{1,ll}.regressor{1,7}.value(8,:));   rr_max1=max(data_in{1,ll}.regressor{1,7}.value(8,:)); % invF1
            rr_min2=min(data_in{1,ll}.regressor{1,8}.value(8,:));   rr_max2=max(data_in{1,ll}.regressor{1,8}.value(8,:)); % invF2
            % range correction
            invF1_corrected=(data_in{1,ll}.regressor{1,7}.value(8,:)-rr_min1)/(rr_max1-rr_min1); %[0,1]
            invF2_corrected=(data_in{1,ll}.regressor{1,8}.value(8,:)-rr_min2)/(rr_max2-rr_min2); %[0,1]
            % mean correction
            invF1_corrected=invF1_corrected-mean(invF1_corrected);
            invF2_corrected=invF2_corrected-mean(invF2_corrected);
            % regressor20: (invF1-invF2)
            data_in{1,ll}.regressor{1,20}.value(8,:)=invF1_corrected-invF2_corrected;
            %             data_in{1,ll}.regressor{1,20}.value(8,zero_pt_ind)=0.0; % fill up '0' for the trivial state (state#1)
            % regressor21: abs(invF1-invF2)
            data_in{1,ll}.regressor{1,21}.value(8,:)=abs(invF1_corrected-invF2_corrected);
            % regressor22- max(invF1-invF2)
            data_in{1,ll}.regressor{1,22}.value=data_in{1,ll}.regressor{1,21}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,22}.value(8,:)=max(invF1_corrected,invF2_corrected);
            % regressor23- conflict btw invF1&invF2
            conflict_sigma=(-1)*log(0.1); % 1% signal remains at the extreme diff. (i.e., max conflict: min conflict = 100:1 )
            data_in{1,ll}.regressor{1,23}.value=data_in{1,ll}.regressor{1,21}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,23}.value(8,:)=exp((-1)*conflict_sigma*((invF1_corrected-invF2_corrected).^2));
            data_in{1,ll}.regressor{1,23}.value(8,zero_pt_ind)=0.0; % fill up '0' for the trivial state (state#1)
            %             % regressor 25 -
            %             mct=data_in{1,ll}.regressor{1,25}.value; % copy the entire structure
            %             % copy the entire structure
            %             data_in{1,ll}.regressor{1,25}.value(7,:) = mean(mct(7,:));
            % regressor26- complexity * max_invM1M2
            data_in{1,ll}.regressor{1,26}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,26}.value(8,:)=data_in{1,ll}.regressor{1,25}.value(7,:).*data_in{1,ll}.regressor{1,22}.value(8,:);
            
            % use complexity (1/2)
            lowcomp_index = find(data_in{1,ll}.regressor{1,25}.value(7,:) == 1);
            highcomp_index = find(data_in{1,ll}.regressor{1,25}.value(7,:) == 2);
            
            %             data_in{1,ll}.regressor{1,27}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            %             data_in{1,ll}.regressor{1,27}.value(7:1:9,:) = zeros(3, size(data_in{1,ll}.regressor{1,22}.value,2));
            %             data_in{1,ll}.regressor{1,27}.value(8,lowcomp_index)=-1.*data_in{1,ll}.regressor{1,22}.value(8,lowcomp_index);
            %
            %             data_in{1,ll}.regressor{1,28}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            %             data_in{1,ll}.regressor{1,28}.value(7:1:9,:) = zeros(3, size(data_in{1,ll}.regressor{1,22}.value,2));
            %             data_in{1,ll}.regressor{1,28}.value(8,highcomp_index)=1.*data_in{1,ll}.regressor{1,22}.value(8,highcomp_index);
            %
            %             data_in{1,ll}.regressor{1,25}.value(8,lowcomp_index)=-1;
            %             data_in{1,ll}.regressor{1,25}.value(8,highcomp_index)=1;
            
            data_in{1,ll}.regressor{1,27}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,27}.value(7:1:9,:) = zeros(3, size(data_in{1,ll}.regressor{1,22}.value,2));
            data_in{1,ll}.regressor{1,27}.value(8,lowcomp_index)=1.*data_in{1,ll}.regressor{1,22}.value(8,lowcomp_index);
            
            data_in{1,ll}.regressor{1,28}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,28}.value(7:1:9,:) = zeros(3, size(data_in{1,ll}.regressor{1,22}.value,2));
            data_in{1,ll}.regressor{1,28}.value(8,highcomp_index)=-1.*data_in{1,ll}.regressor{1,22}.value(8,highcomp_index);
            
            data_in{1,ll}.regressor{1,25}.value(7,lowcomp_index)=1;
            data_in{1,ll}.regressor{1,25}.value(7,highcomp_index)=2;
            
            %             a = length(data_in{1,ll}.regressor{1,25}.value(7,:));
            %             data_in{1,ll}.regressor{1,25}.value(7,1:floor(a/4)) = ones(1, floor(a/4));
            
            data_in{1,ll}.regressor{1,26}.value=data_in{1,ll}.regressor{1,22}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,26}.value(8,lowcomp_index)=1*data_in{1,ll}.regressor{1,22}.value(8,lowcomp_index);
            data_in{1,ll}.regressor{1,26}.value(8,highcomp_index)=2*data_in{1,ll}.regressor{1,22}.value(8,highcomp_index);
            
            data_in{1,ll}.regressor{1,29}.value=data_in{1,ll}.regressor{1,7}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,29}.value(8,lowcomp_index)=1*data_in{1,ll}.regressor{1,7}.value(8,lowcomp_index);
            data_in{1,ll}.regressor{1,29}.value(8,highcomp_index)=2*data_in{1,ll}.regressor{1,7}.value(8,highcomp_index);
            
            data_in{1,ll}.regressor{1,30}.value=data_in{1,ll}.regressor{1,8}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,30}.value(8,lowcomp_index)=1*data_in{1,ll}.regressor{1,8}.value(8,lowcomp_index);
            data_in{1,ll}.regressor{1,30}.value(8,highcomp_index)=2*data_in{1,ll}.regressor{1,8}.value(8,highcomp_index);
            
            behav = [];
            for sei = 1 : length(data_in{1,1}.HIST_block_condition)
                behav = [behav data_in{1,1}.HIST_block_condition{1,sei}(2,:)];
            end
            val_u = [];
            for bei = 1 : length(behav)
                if behav(bei)<3
                    val_u = [val_u -4 -4 -4];
                else
                    val_u = [val_u 4 4 4];
                end
            end
            data_in{1,ll}.regressor{1,31}.value=data_in{1,ll}.regressor{1,25}.value; % copy the entire structure
            data_in{1,ll}.regressor{1,31}.value(7,:)=val_u;
            
            data_in{1,ll}.regressor{1,32}.value=data_in{1,ll}.regressor{1,7}.value;
            data_in{1,ll}.regressor{1,33}.value=data_in{1,ll}.regressor{1,7}.value;
            data_in{1,ll}.regressor{1,34}.value=data_in{1,ll}.regressor{1,7}.value;
            ind_0 = find(data_in{1,ll}.regressor{1,7}.value(8,:) ==0);
            ind_a = 1:1:length(data_in{1,ll}.regressor{1,7}.value(8,:));
            %             data_in{1,ll}.regressor{1,32}.value=data_in{1,ll}.regressor{1,22}.value;
            
            low_tol = quantile (data_in{1,ll}.regressor{1,22}.value(8, setdiff(ind_a, ind_0)), 0.5);
            high_tol = quantile (data_in{1,ll}.regressor{1,22}.value(8, setdiff(ind_a, ind_0)), 0.5);
            data_in{1,ll}.regressor{1,32}.value(8, :) = zeros(1,length(data_in{1,ll}.regressor{1,32}.value(8, :)));
            data_in{1,ll}.regressor{1,32}.value(8, find(data_in{1,ll}.regressor{1,22}.value(8,:) < low_tol)) = -2;
            data_in{1,ll}.regressor{1,32}.value(8, find(data_in{1,ll}.regressor{1,22}.value(8,:) > high_tol)) = 2;
            data_in{1,ll}.regressor{1,32}.value(8, ind_0) = zeros(1,length(ind_0));
            
            ind_0 = find(data_in{1,ll}.regressor{1,7}.value(8,:) ==0);
            low_ind_rel =  find(data_in{1,ll}.regressor{1,32}.value(8,:) == -2);
            tol_ = median(data_in{1,ll}.regressor{1,32}.value(8, setdiff(ind_a,ind_0)));
            low_ind_rel =  setdiff(find(data_in{1,ll}.regressor{1,32}.value(8,:)<tol_), ind_0);
            high_ind_rel = find(data_in{1,ll}.regressor{1,32}.value(8,:) == 2);
            high_ind_rel =  setdiff(find(data_in{1,ll}.regressor{1,32}.value(8,:)>tol_), ind_0);
            low_ind_comp = lowcomp_index;
            high_ind_comp = highcomp_index;
            low_ind_uncer = find(data_in{1,ll}.regressor{1,31}.value(7,:) == -4);
            high_ind_uncer = find(data_in{1,ll}.regressor{1,31}.value(7,:) == 4);
            % REL_COMP
            data_in{1,ll}.regressor{1,33}.value(8, :) = zeros(1,length(data_in{1,ll}.regressor{1,33}.value(8, :)));
            data_in{1,ll}.regressor{1,33}.value(8, intersect(low_ind_comp, low_ind_rel)) = -4;
            data_in{1,ll}.regressor{1,33}.value(8, intersect(low_ind_comp, high_ind_rel)) = -2;
            data_in{1,ll}.regressor{1,33}.value(8, intersect(high_ind_comp,low_ind_rel)) = 2;
            data_in{1,ll}.regressor{1,33}.value(8, intersect(high_ind_comp,high_ind_rel)) = 4;
            
            
            % REL_UNCER
            data_in{1,ll}.regressor{1,34}.value(8, :) = zeros(1,length(data_in{1,ll}.regressor{1,34}.value(8, :)));
            data_in{1,ll}.regressor{1,34}.value(8, intersect(low_ind_comp, low_ind_uncer)) = -4;
            data_in{1,ll}.regressor{1,34}.value(8, intersect(low_ind_comp, high_ind_uncer)) = -2;
            data_in{1,ll}.regressor{1,34}.value(8, intersect(high_ind_comp,low_ind_uncer)) = 2;
            data_in{1,ll}.regressor{1,34}.value(8, intersect(high_ind_comp,high_ind_uncer)) = 4;
            
            
        end
    end
    
    
end


%% fitness values (to be maximized)
fitness_val(pop_id)=1e6*1.0/Sum_NegLogLik;
%     fitness_val(pop_id)=(-1.0)*Sum_NegLogLik;
Sum_NegLogLik_val(pop_id)=Sum_NegLogLik;
%     disp(sprintf('    = fitness value: %04.4f (NegLogLik: %f)',max(fitness_val),min(Sum_NegLogLik)));


%% function returns
if(mode.out==0)
    data_out=OBS;
end
if(mode.out==1)
    data_out=Sum_NegLogLik_val(pop_id,1);
end
if(mode.out==99) % return SBJ structure added with regressors
    data_out=data_in;
end


% function ends
end






