function [out]=batch_model_arbitration_regressor_gen_v5(job_opt)

% clear all
% close all

% plug-in
warning('off')
list_sbj_included=job_opt.list_sbj_included;
%


path0='/home/swlee/fmri_arbitration/';
seed_path_result=[path0 'result_save/'];
save_path_result=[path0 '\result_simul\'];
save_for_SPM=['C:\DATA_fmri\uncertainty_arbitration\regressors_contrasts\'];
save_path_neuroecon=[path0 'regressors_contrasts/'];



% 1. Behavioral data
% LIST_SBJ={'david', 'DeDe', 'rosemary', 'Boyu', 'melissa', 'Rehevolew', 'joel', 'clarke', 'angela', 'william', 'josephine'}; % (good in pre which is mostly habitual - rosemary, melissa)

% 2. behavioral + fmri data (for map config, see SIMUL_arbitraion_fmri2.m)
% [note] 'Oliver' uses an old map. the rest of them use a new map.
LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
LIST_sbj_map_type=[1*ones(1,12) 2*ones(1,50)]; %1:'sangwan2012b', 2:'sangwan2012c'

% regressor list
% [CAUTION] DO NOT change the order!!!
% [NOTE] if "TYPE_REGRESSOR" changed, change "param_regressor_type_cue_abs_pos_in_design_mat" accordingly!!!
LIST_REGRESSOR={'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2', 'Qfwd', 'Qsarsa', 'Qarb', 'dQbwdEnergy', 'dQbwdMean', 'duncertaintyM1', 'dinvFanoM1'};
TYPE_REGRESSOR=[1 1, 1.5 1.5, 1.5 1.5, 1.5 1.5, 2 2, 1.5 1.5 1.5 1.5 1.5 1.5 1.5]; % 1: parametric modulation (0-duration), 1.5:parmetric modulation (non-zero duration), 2: extra continuous parametric modulation (TR-fixed) - this will be used by "dummy" regressor.
row_mat=[7 7 8 8 8 8 8 8 7 7 7 7 7 7 7 8 8]; % from which row in the SBJ{}.regressor matrix the signal needs to be extracted. e.g., uncertainty of 0 prediction error




%% OPTION - subject 
% [note] DO NOT USE sbj#[20] - he pressed wrong buttons in session1,2, so need to shrink all the SBJ matrix size by deleting the session#1,2
list_sbj_included=[2:1:19 21:1:24]; 

%% OPTION - model optimization
option_optimizing_model=0; % 0: optimizing the model for each sbj, 1: for all sbj, 2: do not optimize; load saved model
mode.USE_FWDSARSA_ONLY=0; % 0: arbitration, 1: use fwd only, 2: use sarsa only
mode.USE_BWDupdate_of_FWDmodel=1; % 1: use the backward update for goal-directed model (fwd model), 0: do not use
mode.DEBUG_Q_VALUE_CHG=0; % Debug option 1: show Q-value before/after whenever there is a goal change.
mode.path_ext=path0;
mode.total_simul=20; % # of total simulation repetition per subject
mode.simul_process_display=0; % 1: display model's process, 0: no diplay
mode.experience_sbj_events=[0 1]; % [pre main] - 1: experience exactly the same events(decision,state) as subjects. 0: model's own experience
mode.max_iter=2; % maximum iteration for optimization
% mode.out=1; % 1: normal evaluation mode, 99: regressor added to the SBJ, 0: debug mode






%% OPTION - Regressor arrangement
% 'SPE', 'RPE', 'uncertaintyM1', 'uncertaintyM2', 'meanM1', 'meanM2', 'invFanoM1', 'invFanoM2', 'weigtM1', 'weigtM2'
% should add the regressors in the order of importance
param_regressor_type_cue={'uncertaintyM1', 'invFanoM1', 'uncertaintyM2', 'invFanoM2', 'duncertaintyM1', 'dinvFanoM1', 'dQbwdEnergy', 'Qarb', 'SPE', 'RPE', 'weigtM1'}; 
reg_type_go_first=[1.5 1 2]; % [CAUTION] The order should match with 'param_regressor_type_cue'.   [CAUTION] type"2" should go always last!!!
Is_save_files_local=0; % save optimization parameters and regressor files
Is_save_files_cluster=1; % save optimization parameters and regressor files








%% initialization
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
% make a regressor index for 2nd parametric modulations (dummy)
ind_regressor_type_dummy.ind_reg=ind_regressor_total(find(type_regressor==2));
reg_cnt=reg_cnt+length(ind_regressor_type_dummy.ind_reg);   param_regressor_type_cue_abs_pos_in_design_mat=[param_regressor_type_cue_abs_pos_in_design_mat [param_regressor_type_cue_abs_pos_in_design_mat(end)+1:1:reg_cnt]];
for j=1:1:length(ind_regressor_type_dummy.ind_reg)
    ind_regressor_type_dummy.name{1,j}=LIST_REGRESSOR{1,ind_regressor_type_dummy.ind_reg(j)};
end
% (2) 
ind_regressor_type_base{1,1}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(1)));
ind_regressor_type_base{1,2}.abs_pos_in_design_mat=param_regressor_type_cue_abs_pos_in_design_mat(find(type_regressor==reg_type_go_first(2)));
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
    LIST_SBJ_included{1,k}=LIST_SBJ{1,ind_sbj_included(k)};
end
for i=1:1:num_sbj_included %=1. process only 1 subject
    
    SBJ{1,i}.name=LIST_SBJ{1,ind_sbj_included(i)};
    
    % 'pre' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,ind_sbj_included(i)} '_pre_info.mat'];
    file_name_full=[mode.path_ext '\result_save\' file_name];
    load(file_name_full);
    SBJ{1,i}.HIST_block_condition_pre=HIST_block_condition;
    SBJ{1,i}.HIST_behavior_info_pre=HIST_behavior_info;
    
    % 'fmri' file load: HIST_block_condition{1,session_ind}, HIST_behavior_info{1,session_ind}
    file_name=[LIST_SBJ{1,ind_sbj_included(i)} '_fmri_info.mat'];
    file_name_full=[mode.path_ext '\result_save\' file_name];
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

param_init=[0.3, 15.0, 4, 2.67, 11.7, 10.0, 0.15, 0.1];
param_BoundL=[0.2, 10, 0.5*param_init(3:1:6), 0.01, 0.01]; %[0.01, 3, 1.0, 2.0, 0.5, 5, 0.01, 0.01];
param_BoundU=[1.0/2, 40/2, 2*param_init(3:1:6), 0.5, 0.2]; %[1, 20, 6, 10, 4, 30, 0.5, 0.2];

% ## (way1-each) optimizing for *each* subject and plug the result into each SBJ structure
if(option_optimizing_model==0)    
    for ind_sbj=1:1:size(SBJ,2)
        clear SBJ_test;
        SBJ_test{1,1}=SBJ{1,ind_sbj};        
        disp('############################################')
        disp(['#### optimizing RL-arbitrator for ' sprintf('SBJ#%02d...',ind_sbj)]);
        disp('############################################')
        % [1] model optimization
        mode.out=1;
        myFunc_bu = @(x) eval_ArbitrationRL2(x, SBJ_test, mode); % define a new anonymous function
        [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter));   % X0,LB,UB
        % [2-1] add regressor vector to SBJ
        mode.out=99;
        SBJ_test=eval_ArbitrationRL2(model_BayesArb.param,SBJ_test,mode);
        % [3] Save
        model_BayesArb.mode=mode;
        SBJ_test{1,1}.model_BayesArb=model_BayesArb;
        SBJ{1,ind_sbj}=SBJ_test{1,1};        
        save_file_name=['SBJ_structure.mat'];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' SBJ'])
        end
    end
    option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
end

if(option_optimizing_model==1)
    % ## (way2-batch) optimizing for *all* subjects and plug the result into each SBJ structure
    % [1] model optimization
    mode.out=1;
    myFunc_bu = @(x) eval_ArbitrationRL2(x, SBJ, mode); % define a new anonymous function
    [model_BayesArb.param, model_BayesArb.val]=fminsearchbnd(myFunc_bu, param_init, param_BoundL, param_BoundU, optimset('Display','iter','MaxIter',mode.max_iter));   % X0,LB,UB
    % [2-1] add regressor vector to SBJ
    mode.out=99;
    SBJ=eval_ArbitrationRL2(model_BayesArb.param,SBJ,mode);
    % [3] save
    model_BayesArb.mode=mode;
    for ind_sbj=1:1:size(SBJ,2) % plug in a identical parameter (because this is batch)
        SBJ{1,ind_sbj}.model_BayesArb=model_BayesArb;
    end
    save_file_name=['SBJ_structure.mat'];
    if(Is_save_files_local==1)
        eval(['save ' save_path_result save_file_name ' SBJ'])
    end
    option_optimizing_model=2; % and then write regressors to SBJ structure based on this optimized parameter
end

if(option_optimizing_model==2)
    load_file_name=['SBJ_structure.mat'];
    eval(['load ' save_path_result load_file_name])
    % regressor part deleting and regenerating.
    for ff=1:1:length(list_sbj_included)
        disp(sprintf('- writing regressor to SBJ structure (SBJ%02d)...',list_sbj_included(ff)));
        SBJ0{1,1}=SBJ{1,ff};                
        if(isfield(SBJ0{1,1}, 'regressor')==1)
            SBJ0{1,1}=rmfield(SBJ0{1,1},'regressor'); %remove the regressor field
        end
        mode.out=99;
        model_BayesArb.param=SBJ0{1,1}.model_BayesArb.param;
        SBJ0=eval_ArbitrationRL2(model_BayesArb.param,SBJ0,mode); % refresh and add the regressor part
        SBJ1{1,ff}=SBJ0{1,1};
    end
    SBJ=SBJ1;
end



%test mode
if(option_optimizing_model==3)
    model_BayesArb.param=param_init;
    mode.out=99;
    SBJ=eval_ArbitrationRL2(model_BayesArb.param,SBJ,mode);
end




%% Create regressors
% state. 0.5: fixation mark on, 1: S1, 2: S2, 3: S3, 4: S4, 5: S5,
% 6(+/-)0.1: O1(with win/lost msg), 7(+/-)0.1: O2(with win/lost msg), 8(+/-)0.1: O3(with win/lost msg), 9: O4,
% 10:A1, 11:A2, 20: a short blank page display, -99:fail to choose in time limit, (-) when display off

for jj2=1:1:size(SBJ,2)        % each subject
    disp(sprintf('##### creating regressor structures (sbj%02d/%02d) #######',jj2,size(SBJ,2)));
    for kk2=1:1:size(SBJ{1,jj2}.HIST_behavior_info,2)  % each main session
        
        
        
        mat_work=SBJ{1,jj2}.HIST_event_info{1,kk2};
        num_tot_events=size(SBJ{1,jj2}.HIST_behavior_info{1,kk2},1);
              
        
        
        %% 1. Regressor cue presentation - with parametric modulation (timing: stimulus onset)
        % [0-duration] SPE, RPE
        % [RT-duration] Q_fwd, Q_sarsa
        
        % regressor generation for each main session and save it to a single file that is compatible with SPM
        ind_reg=0;  % corresponds to the size of the structure
        ind_reg_abs=0; % actual number of regressors (including parametric)
        durations={};
        onsets={};
        names={};
        pmod=struct('name',{},'param',{},'poly',{});
        
        
        use_model_regressor_cue=1;
        
        % (1) durations, name, onset        
        [tmp col_on]=find((mat_work(7,:)==1)|(mat_work(7,:)==2)|(mat_work(7,:)==3)|(mat_work(7,:)==4)|(mat_work(7,:)==5)...
            |(mat_work(7,:)==5.9)|(mat_work(7,:)==6.1)|(mat_work(7,:)==6.9)|(mat_work(7,:)==7.1)|...
            (mat_work(7,:)==7.9)|(mat_work(7,:)==8.1)|(mat_work(7,:)==9));
        if(length(col_on)~=num_tot_events*3)
            error('-ERROR: variable ''mat_work'' missed some event extraction. check!')
        end
        
        RT_mat=[];  onset_mat=[];
        prev_trial=0;  show_n_th_times_t=0;        
        param_mat{1,1}=zeros(length(ind_regressor_type_base{1,1}.ind_reg),length(col_on));
        param_mat{1,2}=zeros(length(ind_regressor_type_base{1,2}.ind_reg),length(col_on));
        
        for ll2=1:1:length(col_on)
            
            if(ll2<length(col_on)) % usual case
                pt_on=mat_work(4,col_on(ll2));
                col_off=col_on(ll2)+1;
%                 col_off=col_on(ll2)-1+find(mat_work(7,[col_on(ll2):1:(col_on(ll2)+2)])==0.5); % find the next fixation mark presentation
                pt_off=mat_work(4,col_off);
                RT=pt_off-pt_on;
            else % last event in the session is the outcome presentation
                RT=2.0;
            end
            RT_mat=[RT_mat RT];
            onset_t=mat_work(4,col_on(ll2));
            onset_mat=[onset_mat onset_t];
            
            % fill out regresssor values
            if(use_model_regressor_cue==1)
                
                % regressor type1 
                for nn=1:1:length(ind_regressor_type_base{1,1}.ind_reg)
                    mysession=kk2;
                    myblock=mat_work(1,col_on(ll2)); % block in session
                    mytrial=mat_work(2,col_on(ll2)); % trial in block
                    mytrial_s=mat_work(3,col_on(ll2)); % trial_s in trial (arbitration index: 1 at the second stage, 2 at the third stage)
                    mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(nn)}.value(1:4,:);
                    identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                    col_event=find(identity_tmp==0);
                    param_mat{1,1}(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(nn)}.value(row_mat(ind_regressor_type_base{1,1}.ind_reg(nn)),col_event);
                    
                end
                
                
                % regressor type2
                for nn=1:1:length(ind_regressor_type_base{1,2}.ind_reg)
                    mysession=kk2;
                    myblock=mat_work(1,col_on(ll2)); % block in session
                    mytrial=mat_work(2,col_on(ll2)); % trial in block
                    mytrial_s=mat_work(3,col_on(ll2)); % trial_s in trial
                    mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(nn)}.value(1:4,:);
                    identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                    col_event=find(identity_tmp==0);
                    param_mat{1,2}(nn,ll2)=SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(nn)}.value(row_mat(ind_regressor_type_base{1,2}.ind_reg(nn)),col_event);
                end
                
            end
            
        end        
        
        
        % 1st regressors
        ind_reg=ind_reg+1;        
        onsets{1,ind_reg}=onset_mat;
        if(reg_type_go_first(1)==1)            names{1,ind_reg}=['Cue_{0T}'];  durations{1,ind_reg}=0;        end
        if(reg_type_go_first(1)==1.5)        names{1,ind_reg}=['Cue_{1T}'];    durations{1,ind_reg}=RT_mat;       end
        
        ind_reg_abs=1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors
        
        
        % (2) pmod: how many times each cue presented
        if(use_model_regressor_cue==1)
            for nn=1:1:length(ind_regressor_type_base{1,1}.ind_reg)
                pmod(1,ind_reg).name{1,nn}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,1}.ind_reg(nn)}.name];
                pmod(1,ind_reg).poly{1,nn}=1;
                pmod(1,ind_reg).param{1,nn}=param_mat{1,1}(nn,:);
                ind_reg_abs=ind_regressor_type_base{1,1}.abs_pos_in_design_mat(nn);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,nn};                
            end
        end
        
        
        % 2nd regressors
        ind_reg=ind_reg+1;        
        onsets{1,ind_reg}=onset_mat;
        if(reg_type_go_first(2)==1)            names{1,ind_reg}=['Cue_{0T}'];  durations{1,ind_reg}=0;        end
        if(reg_type_go_first(2)==1.5)        names{1,ind_reg}=['Cue_{1T}'];    durations{1,ind_reg}=RT_mat;       end
        ind_reg_abs=ind_regressor_type_base{1,1}.abs_pos_in_design_mat(end)+1;      list_name_for_contrast{1,ind_reg_abs}=names{1,ind_reg}; % add to the global list of regresssors
        
        % (2) pmod: how many times each cue presented
        if(use_model_regressor_cue==1)            
            ind_reg_param=0;
            for nn=1:1:length(ind_regressor_type_base{1,2}.ind_reg)
                ind_reg_param=ind_reg_param+1;
                pmod(1,ind_reg).name{1,ind_reg_param}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_base{1,2}.ind_reg(ind_reg_param)}.name];
                pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                pmod(1,ind_reg).param{1,ind_reg_param}=param_mat{1,2}(ind_reg_param,:);
                ind_reg_abs=ind_regressor_type_base{1,2}.abs_pos_in_design_mat(ind_reg_param);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};                
            end
            % Reserve spaces for extra regressors (will come from the dummy design matrix)
            for nn=1:1:length(ind_regressor_type_dummy.ind_reg)                
                ind_reg_param=ind_reg_param+1;
                pmod(1,ind_reg).name{1,ind_reg_param}=ind_regressor_type_dummy.name{1,nn};
                pmod(1,ind_reg).poly{1,ind_reg_param}=1;
                pmod(1,ind_reg).param{1,ind_reg_param}=rand(1,length(param_mat{1,2}(1,:)));
                ind_reg_abs=ind_regressor_type_dummy.abs_pos_in_design_mat(nn);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,ind_reg_param};                
            end
        end
        
        
        onset_event=onsets;
        
        
        % (3) Saving normal regressor file
        tot_num_myregressor=length(list_name_for_contrast);
        save_file_name=['Regressor--' SBJ{1,jj2}.name '_sess' sprintf('%02d.mat',kk2)];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' durations names onsets pmod'])
            eval(['save ' save_for_SPM save_file_name ' durations names onsets pmod'])
        end
        if(Is_save_files_cluster==1)
            eval(['save ' save_path_neuroecon save_file_name ' durations names onsets pmod'])
        end
        
        
        
        %% 2. Extra Regressors independent of cue presentation (timing: TR)        
        % The regressors will be saved in separate .mat file.
        
        % regressor generation for each main session and save it to a single file that is compatible with SPM
        ind_reg=0;  % corresponds to the size of the structure
        ind_reg_abs=0; % actual number of regressors (including parametric)
        durations={};
        onsets={};
        names={};
        pmod=struct('name',{},'param',{},'poly',{});
                
        
        ind_reg=ind_reg+1;
        TR_CBIC=2.78; % (sec)                        
        % (1) duration,onset,name
        durations{1,ind_reg}=TR_CBIC;
        onsets{1,ind_reg}=[0:TR_CBIC:(mat_work(4,end)+20)];
        names{1,ind_reg}=['Dummy'];
        % the name will NOT be added to "list_name_for_contrast" because it is a dummy regressor.                
        length_reg=length(onsets{1,1});
        % (2-1) determine the regressor values
        regressor_extra=zeros(length(ind_regressor_type_dummy.ind_reg),length_reg);
        for nn=1:1:length(ind_regressor_type_dummy.ind_reg) % for each regressor            
            % find the previous onset point closest to the current time            
            t_scan=0;
            for i_scan=1:1:length_reg
                % find the recent onset (event) time on which the regressor value has been updated.                
                test_mat=onset_event{1,1}-t_scan;
                test_ii=find(test_mat<0);
                if(length(test_ii)==0) % first a few scans before the event
                    myscan_time=onset_event{1,1}(1);
                else % use the value of previous event
                    myscan_time=onset_event{1,1}(min(length(onset_event{1,1}),test_ii(end)));
                end
                % find the corresponding {block#, trial#, and trial_s#} in event matrix
                mycol=find(abs(mat_work(4,:)-myscan_time)<0.001);
                mysession=kk2; % session
                myblock=mat_work(1,mycol); % block in session
                mytrial=mat_work(2,mycol); % trial in block
                mytrial_s=mat_work(3,mycol); % trial_s in trial (mytrial_s=max(1,mat_work(3,mycol)-1);)                
                % find the corresponding regressor value in a regressor matrix
                if(mytrial_s~=1)
                    mat_work_reg=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(1:4,:);
                    identity_tmp=sum(abs(mat_work_reg-repmat([mysession myblock mytrial mytrial_s]',1,size(mat_work_reg,2))));
                    col_event=find(identity_tmp==0);
                    regressor_extra(nn,i_scan)=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(row_mat(ind_regressor_type_dummy.ind_reg(nn)),col_event);
                else % mytrial_s=1: the first state
                    if(i_scan==1) % for the very first scan, we simply read the first regressor value.
                        regressor_extra(nn,i_scan)=SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.value(row_mat(ind_regressor_type_dummy.ind_reg(nn)),1);
                    else % for the first state in every trial, we simply take the regressor value at t-1 (because there is no update)
                        regressor_extra(nn,i_scan)=regressor_extra(nn,i_scan-1);
                    end                    
                end
                % compute the next scan time
                t_scan=t_scan+TR_CBIC;
            end
        end
        % (2-2) pmod: parametric modulators
        for nn=1:1:length(ind_regressor_type_dummy.ind_reg)
            pmod(1,ind_reg).name{1,nn}=[SBJ{1,jj2}.regressor{1,ind_regressor_type_dummy.ind_reg(nn)}.name];
            pmod(1,ind_reg).poly{1,nn}=1;
            pmod(1,ind_reg).param{1,nn}=regressor_extra(nn,:);
            ind_reg_abs=ind_regressor_type_dummy.abs_pos_in_design_mat(nn);      list_name_for_contrast{1,ind_reg_abs}=pmod(1,ind_reg).name{1,nn};                
        end
        
        % (3) Saving dummy regressor file
        tot_num_myregressor=length(list_name_for_contrast);
        save_file_name=['Regressor_dummy--' SBJ{1,jj2}.name '_sess' sprintf('%02d.mat',kk2)];
        if(Is_save_files_local==1)
            eval(['save ' save_path_result save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
            eval(['save ' save_for_SPM save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
        end
        if(Is_save_files_cluster==1)
            eval(['save ' save_path_neuroecon save_file_name ' durations names onsets pmod ind_regressor_type_dummy'])
        end
        
        
    end
end






%% Saving Contrast file
% [index of my regressors for contrast vector] : total main regressor=6, total regressors=7
clear contrast_spm

total_number_regressor=tot_num_myregressor+6; % # + 6 movements
ind_contrast_vec=0;


% individual : (ex) [0 1 0 0 0 0 0 0]
for ii=1:1:tot_num_myregressor 
    ind_contrast_vec=ind_contrast_vec+1;
    contrast=zeros(1,tot_num_myregressor);
    contrast(1,ii)=1;
    contrast_spm{1,ind_contrast_vec}.name=list_name_for_contrast{1,ii};
    contrast_spm{1,ind_contrast_vec}.vec=contrast;    
end

% % difference : (ex) [0 0 0 0 1 -1 0 0]
% combination_mat=combnk(param_regressor_type_cue_abs_pos_in_design_mat,2);
% for kk=1:1:size(combination_mat,1)
%     % A-B
%     ind_contrast_vec=ind_contrast_vec+1;
%     contrast=zeros(1,total_number_regressor);
%     contrast(1,combination_mat(kk,1))=1;    contrast(1,combination_mat(kk,2))=-1;
%     contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,combination_mat(kk,1)} '>' list_name_for_contrast{1,combination_mat(kk,2)}];
%     contrast_spm{1,ind_contrast_vec}.vec=contrast;
%     % B-A
%     ind_contrast_vec=ind_contrast_vec+1;
%     contrast=zeros(1,total_number_regressor);
%     contrast(1,combination_mat(kk,2))=1;    contrast(1,combination_mat(kk,1))=-1;
%     contrast_spm{1,ind_contrast_vec}.name=[list_name_for_contrast{1,combination_mat(kk,1)} '<' list_name_for_contrast{1,combination_mat(kk,2)}];
%     contrast_spm{1,ind_contrast_vec}.vec=contrast;
% end

if(Is_save_files_local==1)
    eval(['save ' save_path_result 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat'])
    eval(['save ' save_for_SPM 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat'])
end
if(Is_save_files_cluster==1)
    eval(['save ' save_path_neuroecon 'contrast_spm.mat' ' contrast_spm ind_regressor_type_base ind_regressor_type_dummy param_regressor_type_cue_abs_pos_in_design_mat'])
end









%% measure the degree of habit in habitual conditions
% block condition - 1: G(with low uncertainty), 2: G''(with high uncertainty), 3:H(with high uncertainty), 4:H''(with low uncertainty)';
% ### [note]: using "state_action_vec_ref" might be stupid idea. using the
% actual action taken by the model would make more sense!!!
if(1)
    
    state_action_vec_ref=[1 2; 2 1; 3 2; 4 2; 5 1]; % col1:state, col2:corresponding action
    for i=1:1:num_sbj_included
        num_tot_sess=size(SBJ{1,i}.HIST_behavior_info,2);
        sum_mat_percentage=zeros(4,2);
        for i_sess=1:1:num_tot_sess
            condi_to_check=[1 2 3 4]; %for all conditions
            mat_percentage=[];
            for kk=1:1:length(condi_to_check)
                % (1) haibual condition
                row_condi=find(SBJ{1,i}.HIST_behavior_info{1,i_sess}(:,3)==condi_to_check(kk));
                num_total_trial=size(SBJ{1,i}.HIST_behavior_info{1,i_sess},1);
                val_score=0;
                for j=1:1:length(row_condi)
                    state_vec=SBJ{1,i}.HIST_behavior_info{1,i_sess}(row_condi(j),[4:5]);
                    state_action_vec0=state_action_vec_ref(state_vec,2)'; % strong habitual action
                    state_action_vec1=SBJ{1,i}.HIST_behavior_info{1,i_sess}(row_condi(j),[7:8]); % subject's action
                    eval_vec=abs(state_action_vec0-state_action_vec1); % all zero = same actions
                    val_score=val_score+length(find(eval_vec==0));
                end
                mat_percentage=[mat_percentage; [condi_to_check(kk) 100*val_score/(2*length(row_condi))]];
            end
            SBJ{1,i}.HIST_block_condition_habit_score{1,i_sess}=mat_percentage;
            SBJ{1,i}.HIST_block_condition_habit_score_Tag{1,i_sess}='col1: block condition, col2: percentage of habitual action';
            sum_mat_percentage=sum_mat_percentage+mat_percentage;
        end
        SBJ{1,i}.HIST_block_condition_habit_score_mean=sum_mat_percentage/num_tot_sess;
        
    end
    
end






disp('- all done.')






















end