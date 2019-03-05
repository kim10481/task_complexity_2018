function [out]=batch_model_cog_regressor_gen_v5_indi_complex(sis,ori,rvs,m_no) %ori = 0/1; m_no= 1:7; rvs=0/1
% 2017/12 version, 10*6 = 60
% 10 models : CogMod = [0 1 1.12 1.21 2 2.12 2.21 3 3.12 3.21];
% 6 variations : reversedOrNot (2) * 3MFor3QorOriginal (3)

%% JOB SUBMISSION
job_opt.is_on_cluster=0;

addpath(pwd);
addpath([pwd '/modelRLsource']);

job_opt.LIST_SBJ={'Rosebrough', 'Brian', 'Monica', 'Jaymie', 'Arian',...
    'Samuel', 'Christopher', 'Teagan', 'William', 'Conder',...
    'Rebecca', 'Daniel', 'Herman', 'Carrie', 'Abigail',...
    'Danielle', 'Keegan', 'Valerie', 'Devashish', 'Conway',...
    'Jessica', 'Bapat', 'Nat', 'Anand'};
job_opt.LIST_sbj_map_type=[2*ones(1,size(job_opt.LIST_SBJ,2))];

LIST_SBJ={ job_opt.LIST_SBJ{sis:sis} };
sbj_included=[1:1:24];%[1:1:8 10 11 12 14:1:15 17:1:23]; % for 9, 13 and 24, regressor file and func (scans) need to be matched!

if ori == 1
    txt1 = 'ori_';
elseif ori == 2 % 3Q model
    txt1 = ['3Qtauonpsa_'];
    
elseif ori == 3
    txt1 = ['3Qtauonpsa_'];
    
else
    txt1= '3tauonpsa_';
end

if rvs ==0
    
else
    txt1 = ['Reversed_' txt1];
end


CogMod = [0 1 1.12 1.21 2 2.12 2.21];
CogMod = [0 1 1.12 1.21 2 2.12 2.21 3 3.12 3.21];

post_filetext= [txt1 num2str(m_no) '_vMF_Coin']; % postfix of the SBJ_structure file to be saved
job_opt.experience_sbj_events=[1 1]; % [pre main]  +1: experience exactly the same events(decision,state) as subjects. 0: model's own experience -1: use saved setting

job_opt.opt_ArbModel=0; % 0: full arbitrator, 1: invF-based, 2: mean-based,
job_opt.USE_FWDSARSA_ONLY=0; % 0: arbitration, 1: use fwd only, 2: use sarsa only
job_opt.option_optimizing_model=0; % 0: optimizing the model for each sbj, 1: for all sbj, 2: do not optimize; load saved model
job_opt.opt_cogload=CogMod(m_no); % complexity interaction - 0: interaction ignored, 1: interacting with A, 2: interaction with B
job_opt.optimization_max_iter=600; % maximum iteration for optimization
job_opt.num_simul_per_sbj=100; %100 % # of total simulation repetition per subject

%
% [NOTE] need to implement a randomized initialization version (refer to tDCS gen)
%

%% each opt
if(job_opt.option_optimizing_model==0)    % each
    job_opt.post_filetext=['_each_exp_' sprintf('%d%d',job_opt.experience_sbj_events(1),job_opt.experience_sbj_events(2)) post_filetext];
    LIST_ITER=[1:1:length(sbj_included)];
end
if(job_opt.option_optimizing_model==1)    % batch
    job_opt.post_filetext=['_each_exp_' sprintf('%d%d',job_opt.experience_sbj_events(1),job_opt.experience_sbj_events(2)) post_filetext];
    LIST_ITER=[1:1:1];
end

%
% i_job_sub0=1;

%for jj=LIST_ITER
job_opt.list_sbj_included=sbj_included(sis);

warning('off')
list_sbj_included=job_opt.list_sbj_included;
%
if    job_opt.is_on_cluster== 1
    path1='/home/swlee/fmri_cog/';
    path0='/home/swlee/fmri_cog/modelRLsource/';
else
    path1=[pwd '/'];
    path0=[pwd '/modelRLsource/'];
end

seed_path_result=[path0 'result_save/'];
save_path_result=[path0 'result_simul/'];
save_for_SPM=[path1 'regressors_contrasts/'];
save_path_neuroecon=[path1 'for_spm/regressors_contrasts/'];
is_on_cluster=job_opt.is_on_cluster;


% 1. Behavioral data
% LIST_SBJ={'david', 'DeDe', 'rosemary', 'Boyu', 'melissa', 'Rehevolew', 'joel', 'clarke', 'angela', 'william', 'josephine'}; % (good in pre which is mostly habitual - rosemary, melissa)
% mode.map_type=?;

% 2. behavioral + fmri data
% LIST_SBJ=job_opt.LIST_SBJ; %defined above

LIST_sbj_map_type=job_opt.LIST_sbj_map_type; %1:'sangwan2012b', 2:'sangwan2012c'

% regressor list
% [CAUTION] DO NOT change the order!!!
% [NOTE] if "TYPE_REGRESSOR" changed, change "param_regressor_type_cue_abs_pos_in_design_mat" accordingly!!!
LIST_REGRESSOR={'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2', 'Qfwd', 'Qsarsa', 'Qarb', 'dQbwdEnergy', 'dQbwdMean', 'duncertaintyM1', 'dinvFanoM1'};
TYPE_REGRESSOR=[1 1, 1.5 1.5, 1.5 1.5, 1.5 1.5, 2 2, 1.5 1.5 1.5 1.5 1.5 1.5 1.5]; % 1: parametric modulation (0-duration), 1.5:parmetric modulation (non-zero duration), 1.7:parametric modulation (with decision onset)  2: extra continuous parametric modulation (TR-fixed) - this will be used by "dummy" regressor.
row_mat=[7 7 8 8 8 8 8 8 7 7 7 7 7 7 7 8 8]; % from which row in the SBJ{}.regressor matrix the signal needs to be extracted. e.g., uncertainty of 0 prediction error




%% OPTION - subject
% [note] DO NOT USE sbj#[20] - he pressed wrong buttons in session1,2, so need to shrink all the SBJ matrix size by deleting the session#1,2
list_sbj_included= job_opt.list_sbj_included;%[2:1:19 21:1:24];

%% OPTION - model optimization
option_optimizing_model=job_opt.option_optimizing_model; % 0: optimizing the model for each sbj, 1: for all sbj, 2: do not optimize; load saved model
update_SBJ_structure=0; % 0: no update/just read and use, 1: update the changes to the saved SBJ file
mode.opt_ArbModel=job_opt.opt_ArbModel; % 0: full arbitrator, 1: invF-based, 2: mean-based, 3: uncertainty-based arbitrator
mode.USE_FWDSARSA_ONLY=job_opt.USE_FWDSARSA_ONLY; % 0: arbitration, 1: use fwd only, 2: use sarsa only
mode.USE_BWDupdate_of_FWDmodel=1; % 1: use the backward update for goal-directed model (fwd model), 0: do not use
mode.DEBUG_Q_VALUE_CHG=0; % Debug option 1: show Q-value before/after whenever there is a goal change.
mode.path_ext=path0;
mode.total_simul=job_opt.num_simul_per_sbj; % # of total simulation repetition per subject
mode.simul_process_display=0; % 1: display model's process, 0: no diplay
mode.experience_sbj_events=job_opt.experience_sbj_events; % [pre, main]  +1: experience exactly the same events(decision,state) as subjects. 0: model's own experience -1: use saved setting
mode.max_iter=job_opt.optimization_max_iter; % maximum iteration for optimization
% mode.out=1; % 1: normal evaluation mode, 99: regressor added to the SBJ, 0: debug mode

%minho added
mode.opt_cogload = job_opt.opt_cogload;

%% OPTION - Regressor arrangement
% {'SPE', 'RPE', 'uncertaintyM1', 'invFanoM1', 'Qsarsa','Qfwd', 'Qarb','uncertaintyM2', 'invFanoM2', 'duncertaintyM1', 'dinvFanoM1', 'weigtM1'};
% should add the regressors in the order of importance
param_regressor_type_cue={'SPE', 'RPE', 'uncertaintyM1', 'invFanoM1', 'Qsarsa','Qfwd','Qarb','uncertaintyM2', 'invFanoM2', 'weigtM1'};
reg_type_go_first=[1 1.5 2]; % [CAUTION] The order should match with 'param_regressor_type_cue'.   [CAUTION] type"2" should go always last!!!
Do_create_regressors=0;
Is_save_files_local=0; % save optimization parameters and regressor files
Is_save_files_cluster=1; % save optimization parameters and regressor files


%% OPTION - behaviroal analysis & display
Do_behavioral_analysis=[0];
if(Do_behavioral_analysis(1)==1) % dont need to create regressors in behavioral analysis mode!
    Do_create_regressors=0;
end





%% initialization
if(Is_save_files_local==0)    disp('### files will not be saved to your local PC.');     end
if(Is_save_files_cluster==0)    disp('### files will not be saved to the cluster PC.');     end

use_model_regressor_cue=0;  ind_regressor_total=[];   type_regressor=[];    ind_regressor_total_in_design_mat=[];
for ii=1:1:size(param_regressor_type_cue,2) % collect regressor information
    if(strcmp(param_regressor_type_cue{1,ii},'SPE')==1)    use_model_regressor_cue=1;  ind_chk=1;   end
    if(strcmp(param_regressor_type_cue{1,ii},'RPE')==1)    use_model_regressor_cue=1;  ind_chk=2;   end
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM1')==1)    use_model_regressor_cue=1;    ind_chk=3;   end
    if(strcmp(param_regressor_type_cue{1,ii},'uncertaintyM2')==1)    use_model_regressor_cue=1;    ind_chk=4;   end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM1')==1)    use_model_regressor_cue=1;    ind_chk=5;   end
    if(strcmp(param_regressor_type_cue{1,ii},'meanM2')==1)    use_model_regressor_cue=1;    ind_chk=6;   end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM1')==1)    use_model_regressor_cue=1;    ind_chk=7;   end
    if(strcmp(param_regressor_type_cue{1,ii},'invFanoM2')==1)    use_model_regressor_cue=1;    ind_chk=8;   end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM1')==1)    use_model_regressor_cue=1;    ind_chk=9;   end
    if(strcmp(param_regressor_type_cue{1,ii},'weigtM2')==1)    use_model_regressor_cue=1;    ind_chk=10;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qfwd')==1)    use_model_regressor_cue=1;    ind_chk=11;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qsarsa')==1)    use_model_regressor_cue=1;    ind_chk=12;   end
    if(strcmp(param_regressor_type_cue{1,ii},'Qarb')==1)    use_model_regressor_cue=1;    ind_chk=13;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dQbwdEnergy')==1)    use_model_regressor_cue=1;    ind_chk=14;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dQbwdMean')==1)    use_model_regressor_cue=1;    ind_chk=15;   end
    if(strcmp(param_regressor_type_cue{1,ii},'duncertaintyM1')==1)    use_model_regressor_cue=1;    ind_chk=16;   end
    if(strcmp(param_regressor_type_cue{1,ii},'dinvFanoM1')==1)    use_model_regressor_cue=1;    ind_chk=17;   end
    
    % index of regressor in "SBJ" structure
    ind_regressor_total=[ind_regressor_total ind_chk];
    % regressor type
    type_regressor=[type_regressor TYPE_REGRESSOR(ind_chk)];
end
% make a regressor index matrix for 1st parametric modulations (normal)
% (1) make 'param_regressor_type_cue_abs_pos_in_design_mat'
reg_cnt=0; param_regressor_type_cue_abs_pos_in_design_mat=[];
ind_regressor_type_base{1,1}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(1)));
reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,1}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [2:1:reg_cnt]];
ind_regressor_type_base{1,2}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(2)));
reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,2}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+2:1:reg_cnt]];
if(reg_type_go_first(3)==1.7)
    ind_regressor_type_base{1,3}.ind_reg=ind_regressor_total(find(type_regressor==reg_type_go_first(3)));
    reg_cnt=reg_cnt+1+length(ind_regressor_type_base{1,3}.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+2:1:reg_cnt]];
end

% make a regressor index for 2nd parametric modulations (dummy)
ind_regressor_type_dummy.ind_reg=ind_regressor_total(find(type_regressor==2));
reg_cnt=reg_cnt+length(ind_regressor_type_dummy.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+1:1:reg_cnt]];
for j=1:1:length(ind_regressor_type_dummy.ind_reg)
    ind_regressor_type_dummy.name{1,j}=LIST_REGRESSOR{1,ind_regressor_type_dummy.ind_reg(j)};
end
% (2)
ind_regressor_type_base{1,1}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(1)));
ind_regressor_type_base{1,2}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(2)));
if(reg_type_go_first(3)==1.7)
    ind_regressor_type_base{1,3}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(3)));
end
ind_regressor_type_dummy.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==2));


if(sum(abs(param_regressor_type_cue_abs_pos_in_design_mat-sort(param_regressor_type_cue_abs_pos_in_design_mat,'ascend')))~=0)
    error('- ERROR!!!!: the variable ''param_regressor_type_cue_abs_pos_in_design_mat'' should be in ascending order!!!');
end



%% subject data loading

% which subject to be included
% ### READ ONLY ONE SBJ BECAUSE EACH MODEL WILL LEARN EACH SBJ BEHAVIOR.
ind_sbj_included=list_sbj_included;      SUB_ARRAY=list_sbj_included;
num_sbj_included=length(ind_sbj_included);
ind_included=ind_sbj_included;

for k=1:1:num_sbj_included
    LIST_SBJ_included{1,k}=LIST_SBJ{1,num_sbj_included};
end
for i=1:1:num_sbj_included %=1. process only 1 subject
    
    SBJ{1,i}.name=LIST_SBJ{1,num_sbj_included};
    
    % 'pre' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,num_sbj_included} '_pre_info.mat'];
    if(is_on_cluster==0)
        file_name_full=[mode.path_ext 'result_save/' file_name];
    else
        file_name_full=[mode.path_ext 'result_save/' file_name];
    end
    load(file_name_full);
    SBJ{1,i}.HIST_block_condition_pre=HIST_block_condition;
    SBJ{1,i}.HIST_behavior_info_pre=HIST_behavior_info;
    
    % 'fmri' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,num_sbj_included} '_fmri_info.mat'];
    if(is_on_cluster==0)
        file_name_full=[mode.path_ext 'result_save/' file_name];
    else
        file_name_full=[mode.path_ext 'result_save/' file_name];
    end
    load(file_name_full);
    SBJ{1,i}.HIST_behavior_info=HIST_behavior_info;
    SBJ{1,i}.HIST_behavior_info_Tag=HIST_behavior_info_Tag;
    SBJ{1,i}.HIST_event_info=HIST_event_info;
    SBJ{1,i}.HIST_event_info_Tag=HIST_event_info_Tag;
    SBJ{1,i}.HIST_block_condition=HIST_block_condition;
    SBJ{1,i}.HIST_block_condition_Tag=HIST_block_condition_Tag;
    num_tot_session=size(SBJ{1,i}.HIST_behavior_info,2);
    
    SBJ{1,i}.map_type=LIST_sbj_map_type(ind_sbj_included(i));
    
    % [fixing part!!! - for Oliver]
    if(strcmp(SBJ{1,i}.name,'Oliver'))
        for mm=1:1:size(SBJ{1,i}.HIST_event_info,2) % each session
            mat_fixing=SBJ{1,i}.HIST_event_info{1,mm};
            index_delete=zeros(1,size(mat_fixing,2));
            [r_fix, c_fix]= find(mat_fixing(7,:)==9);
            for nn=1:1:length(c_fix)
                % check the previous event
                if(mat_fixing(7, c_fix(nn)-1)~=0.5)
                    index_delete(c_fix(nn))=1;
                end
            end
            [tmp c_keep]=find(index_delete==0);
            mat_fixed=mat_fixing(:,c_keep);
            SBJ{1,i}.HIST_event_info{1,mm}=mat_fixed;
        end
    end
    
    
    % [NOTE] now we have 4 variables: mode.HIST_block_condition_pre, mode.HIST_block_condition, mode.HIST_behavior_info_pre, mode.HIST_behavior_info
    % to read a block condition, use "block_condition=mode.HIST_block_condition{1,session_ind}(2,block_ind); % G:1,G':2,H:3,H':4"
    
    %     swsw_amount_pre = [swsw_amount_pre mode.HIST_behavior_info_pre{1,1}(end,17)];
    tot_amount_earned_main_each_sbj =[];
    for jk=1:1:size(SBJ{1,i}.HIST_behavior_info,2)    tot_amount_earned_main_each_sbj = [tot_amount_earned_main_each_sbj; SBJ{1,i}.HIST_behavior_info{1,jk}(end,17)]; end
    %     swsw_amount_main=[swsw_amount_main tot_amount_earned_main_each_sbj];
end







%% model optimization
% param_in(1): myArbitrator.PE_tolerance_m1 (m1's threshold for zero PE)
% param_in(2): m2_absPEestimate_lr % (before) myArbitrator.PE_tolerance_m2 (m2's threshold for zero PE)
% param_in(3): myArbitrator.A_12
% param_in(x): myArbitrator.B_12 : based on A12
% param_in(4): myArbitrator.A_21
% param_in(x): myArbitrator.B_21 : based on A21
% param_in(5): myArbitrator.tau_softmax/param_sarsa.tau/param_fwd.tau : better to fix at 0,2. This should be determined in a way that maintains softmax values in a reasonable scale. Otherwise, this will drive the fitness value!
% param_in(6): % param_sarsa.alpha/param_fwd.alpha 0.01~0.2 to ensure a good "state_fwd.T" in phase 1


param_init=[0.5   0.2    3    3   3    3    0.2    0.1];%[0.4995   10.7958    2.0531    4.5944   18.1288    5.0023    0.1645    0.0882];
param_BoundL=[0.3   0.05    1    1   1    1    0.1    0.05];% [0.2, 8, 0.8*param_init(3:1:4), 0.12, 0.02];
param_BoundU=[0.75   0.5    10    12   10    12    0.4    0.15];% [0.6, 12, 1.2*param_init(3:1:4), 0.2, 0.14];
mode.boundary_12=0.12;       mode.boundary_21=0.02; % boundary condition : gating fn(1)
switch m_no
    case 1
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
    case 2
    case 3
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
    case 4
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
    case 5
    case 6
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
    case 7
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
    case 8
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
    case 9
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
    case 10
        param_BoundL(4) = param_init(4);
        param_BoundU(4) = param_init(4);
        param_BoundL(6) = param_init(6);
        param_BoundU(6) = param_init(6);
end

mode.param_init=param_init; mode.param_BoundL=param_BoundL; mode.param_BoundU=param_BoundU;
mode.param_length=length(mode.param_init);

% ## (way1-each) optimizing for *each* subject and plug the result into each SBJ structure
if(option_optimizing_model==0)
    for ind_sbj=1:1:size(SBJ,2)
        clear SBJ_test;
        SBJ_test{1,1}=SBJ{1,ind_sbj};
        disp('############################################')
        disp(['#### optimizing RL-arbitrator for ' sprintf('SBJ#%02d...',sis)]);
        disp('############################################')
        % [1] model optimization
        mode.out=1;
        mode.rvs=rvs;
        mode.ori=ori;
        myFunc_bu = @(x) eval_ArbitrationRL5_tauonpsa(x, SBJ_test, mode); % define a new anonymous function: eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        
        [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter));   % X0,LB,UB ('Display','iter','MaxIter',mode.max_iter));
        % [2-1] add regressor vector to SBJ
        mode.out=99;
        SBJ_test=eval_ArbitrationRL5_tauonpsa(model_BayesArb.param,SBJ_test,mode); %: eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        
        % [3] Save
        model_BayesArb.mode=mode;
        SBJ_test{1,1}.model_BayesArb=model_BayesArb;
        SBJ{1,ind_sbj}=SBJ_test{1,1};
        save_file_name=['SBJ_structure' sprintf('_sbj%02d',sis) job_opt.post_filetext '.mat'];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
        if(Is_save_files_cluster==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
    end
    %     option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
end
if(option_optimizing_model==1)
    
    
    SBJrecon = cell(1,length(sbj_included));
    % ## (way2-batch) optimizing for *all* subjects and plug the result into each SBJ structure
    % [0] retrieve intial configuration for skipping pre-training
    %     SBJ_keep=SBJ;
    %     load_file_name=['SBJ_structure(backup,batch,Oct30_4).mat'];
    %     eval(['load ' save_path_result load_file_name]);
    %     for ff1=1:1:length(SBJ_keep)
    %         SBJ_keep{1,ff1}.init_state_fwd=SBJ{1,ff1}.init_state_fwd;    SBJ_keep{1,ff1}.init_state_sarsa=SBJ{1,ff1}.init_state_sarsa;
    %     end
    %     SBJ=SBJ_keep;
    % [1] model optimization
    mode.out = 99;
    for m = 1 : length(sbj_included)
        load_filename = ['SBJ_structure' sprintf('_sbj%02d',sbj_included(m)) job_opt.post_filetext];
        load(['modelRLsource/result_simul/' load_filename]);
        SBJrecon{1,m}=SBJ{1,1};
        
    end
    option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
    SBJ= SBJrecon;
    save_file_name=['SBJ_structure' job_opt.post_filetext '.mat'];
    eval(['save ' save_path_result save_file_name ' SBJ']);
    eval(['save ' save_path_result 'SBJ_structure.mat' ' SBJ']);
end




if(option_optimizing_model==2) % [NOTE] replace initial SBJ = just read SBJ from the "SBJ_structure.mat"
    load_file_name=['SBJ_structure.mat'];
    eval(['load ' save_path_result load_file_name])
    % regressor part deleting and regenerating.
    for ff=1:1:length(sbj_included)
        disp(sprintf('- writing regressor to SBJ structure (SBJ%02d)...',sbj_included(ff)));
        % find my subject in "SBJ" strucure of the 'SBJ_structure.mat' file
        did_find=0;
        for ss=1:1:size(SBJ,2)
            if(strcmp(SBJ{1,ss}.name,job_opt.LIST_SBJ{1,ff})==1)
                SBJ0{1,1}=SBJ{1,ss}; % SBJ : this includes SBJ structure for subjects to be included for this code
                did_find=did_find+1;
            end
        end
        if(did_find~=1)            error('-ERROR:: no correponding subject found in the "SBJ_structure.mat" file!!!');   end
        
        if(isfield(SBJ0{1,1}, 'regressor')==1)
            SBJ0{1,1}=rmfield(SBJ0{1,1},'regressor'); %remove the regressor field
        end
        mode.out=99;
        mode.rvs=rvs;
        mode.ori=ori;
        model_BayesArb.param=SBJ0{1,1}.model_BayesArb.param;
        SBJ0=eval_ArbitrationRL5_tauonpsa(model_BayesArb.param,SBJ0,mode); %: eval_ArbitrationRL2(x, SBJ_test, mode) for full BayesArb
        SBJ1{1,ff}=SBJ0{1,1};
    end
    clear SBJ
    SBJ=SBJ1;
    if(update_SBJ_structure==1)        eval(['save ' save_path_result load_file_name ' SBJ']);  end
end