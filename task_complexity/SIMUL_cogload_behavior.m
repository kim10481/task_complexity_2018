function [output_info]=SIMUL_cogload_behavior(EXP_NAME_IN,index_num,session_opt)
% [NOTE] the goal is changing in each trial.
% (ex) SIMUL_arbitration_fmri('john',1,'pre')
% (ex) SIMUL_arbitration_fmri('john',1,'fmri')
% works for the map type 'sangwan2014b' ONLY.

% testtest
%  EXP_NAME_IN='test2';
%  index_num=1; % : session#
%  session_opt='pre'; %'pre','fmri'

% scanner key device : USB HHSC-2x4-C
% use "key_test.m" to see the actual number of keyboard input.qwe


%% map name
name_map='sangwan2014b';


%% organizing all the cue presentation files before the practice session
% refer to ../result_save/List_subject_info.xlsx and update!
if(strcmp(session_opt,'pre')==1)
    if(index_num==1) % for the first session per each subject
        disp('- now starting to organize the cue presentation files...');
        Info_seed=SIMUL_cogload_fmri_init(name_map);
        disp('- organization of the cue presentation files completed.');
    end
end

%% function option check
okay_to_start=0;
if(strcmp(session_opt,'pre')==1)
    okay_to_start=1;
end
if(strcmp(session_opt,'fmri')==1)
    okay_to_start=1;
end
if(okay_to_start~=1)
    error('- ERROR:: check the "session_opt!"');
end

EXP_NAME=[EXP_NAME_IN '_' session_opt];

% seed_path0=['F:\0-Program\MATLABwork\work\One shot learning']; % desktop
seed_path0=pwd;%['D:\0-program\One shot learning']; % laptop
seed_path=[seed_path0 '\seed\'];
% save_path=[seed_path0 '\results_save\'];

%% check if session file exists, and other conditions
COND_NEW_FILE=1;
file_chk_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_check=[pwd '\result_save\' file_chk_name];
if(exist(file_name_check)==2)
    disp('$$$$$ ERROR: the corresponding file exists. try another name or session number. $$$$$$$$$$');
    COND_NEW_FILE=0;
    return;
end
if(index_num>1)
    file_chk_name2=[EXP_NAME sprintf('_%d.mat',index_num-1)];
    file_name_check2=[pwd '\result_save\' file_chk_name2];
    if(exist(file_name_check2)==0) % if the previous session file does not exist
        disp('$$$$$ ERROR: The previous session file does not exist. Check the previous session number. $$$$$$$$$$');
        COND_NEW_FILE=0;
        return;
    end
end
MAX_SESSION_NUM=6;
if(index_num>MAX_SESSION_NUM) % MAX session number =5 !
    disp(sprintf('$$$$$ ERROR: exceeds the max # of sessions =%d. $$$$$$$$$$',MAX_SESSION_NUM));
    COND_NEW_FILE=0;
    return;
end

warning('off')
close all
output_info=0;



%% options - display

% for test purpose
DO_SHOW_TRANSITION_PROB=0; % display transition probability 

% screen size - fixed
SCREEN_RESOLUTION=[1024 768];

% image size
KEEP_PICTURE_OPT=1; % always 1
DO_TAKE_SNAPSHOT=0;
IMAGE_SIZE=[256 256]; %width,height - for cue presentation
% GOAL_IMG_SIZE=[256 245]; % width,height - for goal state presentation
CHOICESET_IMG_SIZE=[1152 170]; % width,height - for choice availability presentation
COIN_IMG_SIZE=[200 200]; % height in pixel
disp_scale=1.0; % for cue presentation
disp_scale_goalimg=0.5; % for goal image presentation
disp_scale_outcome_msg=0.5; % for outcome msg presentation
disp_scale_scoring=0.7; % for mouse-clicking score submission display
Tot_session=1; % # of total sessions (fixed because this runs for each session)

% text size
text_size_default=20; % font size (don't change)
text_size_rwd_msg=20;
text_size_currency=30;

% background color, [210,210,210,150]; % light gray
BackgroundColor_block_intro=[165,165,165,150]; % gray
BackgroundColor_Cue_page=BackgroundColor_block_intro;
BackgroundColor_Trial_ready_page=BackgroundColor_block_intro;
BackgroundColor_Reward_page=BackgroundColor_block_intro;
COLOR_FIXATION_MARK=[70,70,70,200]; % dark gray


if(strcmp(session_opt,'pre')==1) % pre-session
    Tot_block=8; % # of total blocks (MUST be the multitude of 4(=#of conditions)
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    Tot_block=4; % 16; % # of total blocks (MUST be the multitude of 4(=#of conditions)
end


%% options - display speed, tiral length, ITI, timing,...

% block schedule
        % _Clow_Ulow: 2 actions available:low complexity, state transition prob=(0.9,0.1):low uncertainty 
        % _Chigh_Ulow: 4 actions available:high complexity, state transition prob=(0.9,0.1):low uncertainty 
        % _Clow_Uhigh: 2 actions available:low complexity, state transition prob=(0.5,0.5):high uncertainty 
        % _Chigh_Uhigh: 4 actions available:high complexity, state transition prob=(0.5,0.5):high uncertainty 
range_num_trials_C_U=[[14,16];[14,16];[14,16];[14,16]]; % each row: minmax # of trials of each block type above

time_estimation_trial_sec=16; %(sec)- used to estimate session time when scheduling
time_limit_session_min=17.8; %(min) - rescheduling until the time estimation meets thie criterion (limit)

if(strcmp(session_opt,'pre')==1) % pre-session
    ffw_speed=1.8;%4; % fast forward speed, 1: normal speed
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    ffw_speed=1;%1; % fast forward speed, 1: normal speed
end
% sec_stim_ready=.1; %(sec)
% sec_trial_ready=1.0/ffw_speed; %(sec)
sec_scanner_ready=5/ffw_speed; % sec for scanner stabilization
% sec_block_ready=0.5/ffw_speed; % sec for block ready signal
% sec_stim_display=0.0/ffw_speed;
sec_currency_disp=3/ffw_speed; %(sec)
sec_stim_interval=[1 4]/(ffw_speed);%1.5; %(sec)
sec_trial_interval=[1 3]/(ffw_speed);%1.5; %(sec)
sec_limit_decision=4;%/(ffw_speed); % time limit of decision (sec)
sec_jittered_blank_page=0.15/(ffw_speed); % (sec)
sec_reward_display=2/ffw_speed; % only for 'fmri' option. 1.5sec for 'pre' session



%% key code definition

% fMRI scanner system setting
KEY_L1=49;
KEY_L2=50;
KEY_R2=51;
KEY_R1=52;
KEY_Y=89; %'y'
KEY_N=78; %'n'
KEY_Q=81; %'q'
KEY_T=53; % 't', 5 in desktop, 84 in laptop

% my laptop setting
% KEY_L1=74;
% KEY_L2=75;
% KEY_R2=76;
% KEY_R1=186;
% KEY_Y=89; %'y'
% KEY_N=78; %'n'
% KEY_Q=81; %'q'
% KEY_T=84; % 't', 5 in desktop, 84 in laptop






%%
%% INITIALIZATION STARTS HERE
%%





%% creating the mother map and state
map_opt.transition_prob_seed=[0.9 0.1];
% map_opt.reward_seed=[40 20 10 0];
[myMap N_state N_action N_transition]=Model_Map_Init2(name_map,map_opt);
% create my state
myState=Model_RL_Init(N_state,N_action,N_transition);


%% global variables (for each block)
HIST_event_info0=[];
HIST_behavior_info0=[];
if(index_num==1) % create the image usage matrix if this is the first session

        % _Clow_Ulow: 2 actions available:low complexity, state transition prob=(0.9,0.1):low uncertainty 
        % _Chigh_Ulow: 4 actions available:high complexity, state transition prob=(0.9,0.1):low uncertainty 
        % _Clow_Uhigh: 2 actions available:low complexity, state transition prob=(0.5,0.5):high uncertainty 
        % _Chigh_Uhigh: 4 actions available:high complexity, state transition prob=(0.5,0.5):high uncertainty 

    
    % event   : HIST_event_info{1,session#}
    HIST_event_info_Tag{1,1}='row1 - block#';    HIST_event_info_Tag{2,1}='row2 - trial# (in each block), 0 if outside of the trial';     HIST_event_info_Tag{3,1}='row3 - trial_s#, 0 if outside of the trial_s';
    HIST_event_info_Tag{4,1}='row4 - event time in session';   HIST_event_info_Tag{5,1}='row5 - event time in block';       HIST_event_info_Tag{6,1}='row6 - event time in trial';
    HIST_event_info_Tag{7,1}='row7 - state. 0.5: fixation mark on, 1~11: S1~S11 >> (+/-)0.1: (with win/lost msg), 21-24:L1/R1/L2/R2, 30: a short blank page display, -99:fail to choose in time limit, (-) when display off';
    HIST_event_info_Tag{8,1}='row8 - goal value of the 1st outcome type';    HIST_event_info_Tag{9,1}='row9 - goal value of the 2nd outcome type';    HIST_event_info_Tag{10,1}='row10 - goal value of the 3rd outcome type';

    % behavior : HIST_behavior_info{1,session#}    
    HIST_behavior_info_Tag{1,1}='col1 - block #';    HIST_behavior_info_Tag{1,2}='col2 - trial # (in each block)';
    HIST_behavior_info_Tag{1,3}='col3 - block condition - 1: Clow_Ulow, 2: Chigh_Ulow, 3: Clow_Uhigh, 4: Chigh_Uhigh';
    HIST_behavior_info_Tag{1,4}='col4 - S1';
    HIST_behavior_info_Tag{1,5}='col5 - S2';        
    HIST_behavior_info_Tag{1,6}='col6 - S3';
    HIST_behavior_info_Tag{1,7}='col7 - A1 (action in the first stage)';     
    HIST_behavior_info_Tag{1,8}='col8 - A2 (action in the second stage)'; 
    HIST_behavior_info_Tag{1,9}='col9 - RT(A1)';        
    HIST_behavior_info_Tag{1,10}='col10 - RT(A2)';
    HIST_behavior_info_Tag{1,11}='col11 - onset (S1) from the trial start';      
    HIST_behavior_info_Tag{1,12}='col12 - onset (S2) from the trial start';      
    HIST_behavior_info_Tag{1,13}='col13 - onset (S3) from the trial start';
    HIST_behavior_info_Tag{1,14}='col14 - onset (A1) from the trial start';      
    HIST_behavior_info_Tag{1,15}='col15 - onset (A2) from the trial start';
    HIST_behavior_info_Tag{1,16}='col16 - reward amount at S3';
    HIST_behavior_info_Tag{1,17}='col17 - total amount (total in the current session)';
    HIST_behavior_info_Tag{1,18}='col18 - goal value of the 1st outcome type';
    HIST_behavior_info_Tag{1,19}='col19 - goal value of the 2nd outcome type';
    HIST_behavior_info_Tag{1,20}='col20 - goal value of the 3rd outcome type';
    
else % load image usage matrix to update
    file_imgind_ld_name=[EXP_NAME '_info.mat'];
    file_name_ld=[pwd '\result_save\' file_imgind_ld_name];
    load(file_name_ld);
end



% seed image read
img_set=cell(1,N_state); % including othe output states
img_set_choice_type=cell(1,3);
img_set_reward_type=cell(1,3);
for img_ind=1:1:N_state
    file_full_path=[seed_path sprintf('s%03d.png',img_ind)];
    img_set{1,img_ind}=imresize(imread(file_full_path),[IMAGE_SIZE(2) IMAGE_SIZE(1)]); % get a frame
end
% choice type read (choice availability)
for jj=1:1:2
    file_full_path=[seed_path sprintf('choice_type%d.png',jj)];
    img_set_choice_type{1,jj}=imresize(imread(file_full_path),[CHOICESET_IMG_SIZE(2) CHOICESET_IMG_SIZE(1)]);
end
% reward image (coin) read
for jj=1:1:3
    file_full_path=[seed_path sprintf('o%03d.png',jj)];
    img_set_reward_type{1,jj}=imresize(imread(file_full_path),[COIN_IMG_SIZE(2) COIN_IMG_SIZE(1)]);
end



%% Scheduling (block sequence)
% _Clow_Ulow: 2 actions available:low complexity, state transition prob=(0.9,0.1):low uncertainty
% _Chigh_Ulow: 4 actions available:high complexity, state transition prob=(0.9,0.1):low uncertainty
% _Clow_Uhigh: 2 actions available:low complexity, state transition prob=(0.5,0.5):high uncertainty
% _Chigh_Uhigh: 4 actions available:high complexity, state transition prob=(0.5,0.5):high uncertainty 

criterion=0;
HIST_block_condition_Tag{1,1}='row1: block #';
HIST_block_condition_Tag{2,1}='row2: block condition for each trial - 1: Clow_Ulow, 2: Chigh_Ulow, 3: Clow_Uhigh, 4: Chigh_Uhigh';
HIST_block_condition_Tag{3,1}='row3: reward value for the 1st outcome type';
HIST_block_condition_Tag{4,1}='row4: reward value for the 2nd outcome type';
HIST_block_condition_Tag{5,1}='row5: reward value for the 3rd outcome type';
if(strcmp(session_opt,'pre')==1) % pre-session (5trials for each block)
    ind_block_condi_pre=2; % all actions available with low state trasition uncertainty
    currency_val=round(myMap.max_currency/2); % no currency change, fixed at the mid value
    tmp_block_ind_mat=[1:1:Tot_block]'*ones(1,5);
    HIST_block_trial_index=[reshape(tmp_block_ind_mat',[1 5*Tot_block]);...
        [ind_block_condi_pre*ones(1,5*Tot_block)];...
        [currency_val*ones(3,5*Tot_block)]]; % no currency change
end
if(strcmp(session_opt,'fmri')==1) % fmri-session
    disp('- scheduling start...');
    while(criterion==0)
        HIST_block_trial_index=[]; % (2xtrial#) row1: block#, row2: the block condition of each trial
        block_index=0;
        for block_quad=1:4:Tot_block % four-block scheduling - to minimize overlap
            current_block_seq=randperm(size(range_num_trials_C_U,1)); % G:1, G':2, H:3, H':4
            if(block_quad>1)
                while(current_block_seq(1)==HIST_block_trial_index(end)) % do not allow the same subsequent condition.
                    current_block_seq=randperm(4); % block condition
                end
            end
%             if((block_quad==1)&&(index_num==1)) % the first block of the first session should be G
%                 current_block_seq=[3 3 1 2];
%             end

            for jj=1:1:length(current_block_seq)
                
                block_index=block_index+1;
                tmp=randperm(max(range_num_trials_C_U(jj,:))-min(range_num_trials_C_U(jj,:))+1); % pick integer
                num_trial=tmp(1)+min(range_num_trials_C_U(jj,:))-1;

                %                 condi_habit=(current_block_seq(jj)==3)||(current_block_seq(jj)==4);
                %                 if((rand>0.52)&&(condi_habit)) % add the devaluation index matrix (in the third row)
                %                     deval_mat0=[ones(1,num_trial-3) zeros(1,3)]; % devaluation for the last two trials
                %                 else
                %                     deval_mat0=[ones(1,num_trial)];
                %                 end
                %                 if((block_quad==1)&&(index_num==1)) % no devaluation for the first block of the first session
                %                     deval_mat0=[ones(1,num_trial)];
                %                 end
                %                 HIST_block_trial_index=[HIST_block_trial_index [ones(1,num_trial)*block_index; ones(1,num_trial)*current_block_seq(jj); deval_mat0]];
                HIST_block_trial_index=[HIST_block_trial_index [ones(1,num_trial)*block_index; ones(1,num_trial)*current_block_seq(jj)]];
                
            end
            
        end
        %% add currency to the block schedule
        condi_curr=0;
        while(condi_curr==0)
            mat_currency=floor((myMap.max_currency+1-0.00001)*rand(3,size(HIST_block_trial_index,2)));
            mat_tmp=corr(mat_currency');
            % make sure that outcome values are not correlated (corr<0.05)
            mat_corr=[mat_tmp(1,2) mat_tmp(1,3) mat_tmp(2,3)];
            if(sum(abs(mat_corr)>0.05)==0)
                condi_curr=1;        HIST_block_trial_index=[HIST_block_trial_index; mat_currency];
            end
        end        
        time_est=size(HIST_block_trial_index,2)*time_estimation_trial_sec/60; %min
        criterion=(time_est<time_limit_session_min);
        if(criterion==0)
            disp('- rescheduling to meet the session time limit criterion...')
        end
    end
    
    str=sprintf('- scheduling done. # of trials = %d. Estimated session time = %02.1f. Avg currency corr = %0.3f.  Proceed? (''n'' to quit, proceed otherwise)',size(HIST_block_trial_index,2),time_est, mean(mat_corr)); disp(str);
    WaitSecs(0.5);
    
    [secs, keyCode] = KbPressWait;      [tmp tmp_key_code]=find(keyCode==1);
    if(tmp_key_code==KEY_N) % n pressed
        disp('### Experiment aborted as per user''s request. ###');    return;
    end
    
end
HIST_block_condition{1,index_num}=HIST_block_trial_index;
disp('- proceed to the experiment...');





%%
%% EXPERIMENT STARTS HERE
%%


%% Display initialization
whichScreen = 1;
wPtr  = Screen('OpenWindow',whichScreen);
[screenWidth, screenHeight] = Screen('WindowSize', wPtr);

white = WhiteIndex(wPtr); % pixel value for white
black = BlackIndex(wPtr); % pixel value for black
gray = (white+black)/2;
inc = white-gray;
inc_0=white-black;

imageArray={};


%% starting message
Screen('TextSize',wPtr, text_size_default);
% Screen('TextFont',wPtr, 'Times New Roman');
DrawFormattedText(wPtr, 'Are you ready for the experiment?\n(Press any key to wait for the trigger)', 'center', 'center');
Screen('Flip', wPtr);  
KbWait; % temporarily disabled for test APR 21
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end


%% waiting for the trigger sign from the scanner
if(strcmp(session_opt,'fmri')==1)
    
    DrawFormattedText(wPtr, 'Now waiting for the trigger...', 'center', 'center');
    Screen('Flip',wPtr);
    % Look for trigger pulse
    while 1 % temporarily disabled for test APR 21
        [ keyIsDown, timeSecs, keyCode ] = KbCheck; % or KbCheck([-1]) for checking from all devices
        if keyIsDown
            [tmp tmp_key_code]=find(keyCode==1);
            
            %         if ((tmp_key_code==53)) % trigger
            if ((tmp_key_code==KEY_T)) % trigger
                %         if (KbName(keyCode) == 5)
                break;
            end
            % If the user holds down a key, KbCheck will report multiple events.
            % To condense multiple 'keyDown' events into a single event, we wait until all
            % keys have been released.
            while KbCheck; end
        end
    end
    
end

%% clock-start and then wait for another 5secs until the scanner stabilizes
session_clock_start = GetSecs;
WaitSecs(sec_scanner_ready);


%% block
for block=1:1:Tot_block % each block

    %% I. Stimuli presentation

    zzz=[];
    % block starts
    %     Screen('FillRect',wPtr,BackgroundColor_block_intro);
    %     str_block_intro=sprintf('*** Now block %d starts. ***',block);
    %     DrawFormattedText(wPtr, str_block_intro, 'center', 'center',[0, 0, 0, 255]);
    %     Screen(wPtr, 'Flip');
        block_clock_start = GetSecs;
    %     HIST_event_info0=[HIST_event_info0 [block; 0; 0; (GetSecs-session_clock_start); (GetSecs-block_clock_start); 0.5]]; % event save
    %     WaitSecs(sec_block_ready);

    trial_in_block=0;
    [tmp trial_set]=find(HIST_block_trial_index(1,:)==block);
    for trial=trial_set % each trial

        trial_in_block=trial_in_block+1;
        trial_clock_start = GetSecs;
        
        onset_state_from_trial_start=[];
        onset_action_from_trial_start=[];

        state_sbj=myState;
        map_sbj=myMap;


        %% [1] map preparation

        % block_condition:
        % _Clow_Ulow: 2 actions available:low complexity, state transition prob=(0.9,0.1):low uncertainty
        % _Chigh_Ulow: 4 actions available:high complexity, state transition prob=(0.9,0.1):low uncertainty
        % _Clow_Uhigh: 2 actions available:low complexity, state transition prob=(0.5,0.5):high uncertainty
        % _Chigh_Uhigh: 4 actions available:high complexity, state transition prob=(0.5,0.5):high uncertainty
        
        Tprob_Ulow=[0.9 0.1];   Tprob_Uhigh=[0.5 0.5];
        action_valid=[1 1 1 1]; % L1,R1,L2,R2

        block_condition=HIST_block_trial_index(2,trial);    txt_transition_P='';
        rwd_condition=HIST_block_trial_index(3:5,trial); % 1x3: rwd value for the three types of goal state

        switch block_condition
            case 1
                Tprob0=Tprob_Ulow;
                action_valid=[[1 1 -99 -99]; [1 1 -99 -99]];
                current_choice_type=2;                
            case 2
                Tprob0=Tprob_Ulow;
                action_valid=[[1 1 -99 -99]; [1 1 1 1]];
                current_choice_type=1;
            case 3
                Tprob0=Tprob_Uhigh;
                action_valid=[[1 1 -99 -99]; [1 1 -99 -99]];
                current_choice_type=2;
            case 4
                Tprob0=Tprob_Uhigh;
                action_valid=[[1 1 -99 -99]; [1 1 1 1]];
                current_choice_type=1;
        end
        txt_transition_P=sprintf('P=(%1.1f,%1.1f)',Tprob0(1),Tprob0(2));
        % update T_prob
        map_sbj=Model_Map_update(map_sbj,'T',Tprob0);
        % set Reward
        map_sbj=Model_Map_update(map_sbj,'R',rwd_condition);
        


        while(map_sbj.IsTerminal(state_sbj.state_history(state_sbj.index))==0)

            
            current_state=state_sbj.state_history(state_sbj.index);


            %% [2] display

            % (1) display fixation mark and display off during the jittered interval
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            if(current_state==1)
                %% show currency values in the first stage
                for i_cc=1:1:3                    
                    % coin img
                    sx2=floor(COIN_IMG_SIZE(1)*disp_scale_outcome_msg);       sy2=floor(COIN_IMG_SIZE(2)*disp_scale_outcome_msg);
                    xt = round(screenWidth/2) + sx2*(i_cc-2);    yt = round(screenHeight/2);
                    input_stim2 = Screen('MakeTexture', wPtr, img_set_reward_type{1,i_cc});                    
                    destrect2=[xt-sx2/2,yt-sy2/2,xt+sx2/2,yt+sy2/2];
                    Screen('DrawTexture', wPtr, input_stim2,[],destrect2); 
                    txt_currency=sprintf('%d',rwd_condition(i_cc));
                    % currency text
                    oldTextSize=Screen('TextSize', wPtr ,text_siqze_currency);
                    DrawFormattedText(wPtr, txt_currency, xt-round(text_size_currency/2), destrect2(4)+10, COLOR_FIXATION_MARK); 
                    oldTextSize=Screen('TextSize', wPtr ,oldTextSize);
                end
                % trial start message
                oldTextSize=Screen('TextSize', wPtr ,text_size_currency);
                DrawFormattedText(wPtr, 'Currency in this trial', 'center', destrect2(2)-100, COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
                oldTextSize=Screen('TextSize', wPtr ,oldTextSize);
            else
                DrawFormattedText(wPtr, '.', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
            end
            Screen(wPtr, 'Flip');
            HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                (GetSecs-session_clock_start); (GetSecs-block_clock_start); (GetSecs-trial_clock_start); 0.5; rwd_condition]]; % event save
            if(current_state==1)
                sec_stim_interval0=sec_currency_disp; % longer sec for the currency display
            else
                sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
            end            
            WaitSecs(sec_stim_interval0);
            
            % (2-1) add state image
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            input_stim = Screen('MakeTexture', wPtr, img_set{1,current_state});
            xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
            sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
            destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];            destrect_state=destrect;
            Screen('DrawTexture', wPtr, input_stim,[],destrect);

            % (2-2) add choice set image            
            disp_choice_type=1;
            if(state_sbj.index==1) % 2 choices available for the first stage
                disp_choice_type = 2; 
            end
            if(state_sbj.index==2)
                disp_choice_type = current_choice_type;
            end
            input_stim2 = Screen('MakeTexture', wPtr, img_set_choice_type{1,disp_choice_type});
            sx2=floor(CHOICESET_IMG_SIZE(1)*disp_scale_goalimg);       sy2=floor(CHOICESET_IMG_SIZE(2)*disp_scale_goalimg);
            xpos2=xpos; ypos2=ypos+sy/2+sy2/2+50;
            destrect2=[xpos2-sx2/2,ypos2-sy2/2,xpos2+sx2/2,ypos2+sy2/2];            
            Screen('DrawTexture', wPtr, input_stim2,[],destrect2);
            % (for test only) show state-transition P
            if(DO_SHOW_TRANSITION_PROB==1)
                oldTextSize=Screen('TextSize', wPtr ,30);
                DrawFormattedText(wPtr, txt_transition_P, xpos, 50, COLOR_FIXATION_MARK);
                oldTextSize=Screen('TextSize', wPtr ,oldTextSize);
            end

            
            
%             % testtest for TEST only - APR 27, 2012
%             if(block_condition==1)                txt_test='G  (0.9, 0.1)';            end
%             if(block_condition==2)                txt_test='G'' (0.5, 0.5)';            end
%             if(block_condition==3)                txt_test='H  (0.5, 0.5)';            end
%             if(block_condition==4)                txt_test='H'' (0.9, 0.1)';            end
%             txt_test2=sprintf('current state=%d,',current_state);
%             DrawFormattedText(wPtr, [txt_test, ', ' txt_test2], 100, 200);
            
            % (2-3) display on
            Screen(wPtr, 'Flip'); % display on
            clock_time_limit_start=GetSecs;
            onset_state_from_trial_start=[onset_state_from_trial_start (GetSecs-trial_clock_start)];
            HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_state_from_trial_start(end); current_state; rwd_condition]]; % event save

            %% [3] get chioce and update            
            decision_made=0;
            while(~decision_made)
                [secs, keyCode] = KbPressWait([], clock_time_limit_start+sec_limit_decision); % if no keyboard in time limit, then go ahead. if pressed earlier, then go ahead.
                onset_action_from_trial_start=[onset_action_from_trial_start (GetSecs-trial_clock_start)];
                state_sbj.RT(state_sbj.index)=GetSecs-clock_time_limit_start;
                [tmp tmp_key_code]=find(keyCode==1);
                if(tmp_key_code == KEY_L1*action_valid(state_sbj.index,1))
                    decision_made=1;
                    state_sbj.action_history(state_sbj.index)=1;
                end
                if(tmp_key_code == KEY_R1*action_valid(state_sbj.index,2))
                    decision_made=1;
                    state_sbj.action_history(state_sbj.index)=2;
                end
                if(tmp_key_code == KEY_L2*action_valid(state_sbj.index,3))
                    decision_made=1;
                    state_sbj.action_history(state_sbj.index)=3;
                end
                if(tmp_key_code == KEY_R2*action_valid(state_sbj.index,4))
                    decision_made=1;
                    state_sbj.action_history(state_sbj.index)=4;
                end
                if(tmp_key_code==KEY_Q) % 'q' pressed for aborting
                    state_sbj.action_history(state_sbj.index)=ceil(2*rand); % random select if fail to make a decision
                    HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                        (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); -99; rwd_condition]]; % event save
                    qwe;
                    break;
                    decision_made=1;                    
                end                
                % check the time limit !@#$                
                if((state_sbj.RT(state_sbj.index)>sec_limit_decision)&&(decision_made==0)) % no decision made in time limit
                    state_sbj.action_history(state_sbj.index)=ceil(2*rand); % random select if failed to make a decision in time
                    decision_made=1;            Is_bet=1;
                    HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                        (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); -99; rwd_condition]]; % event save
                else % in time, and decision made
                    if(decision_made==1)
                        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                            (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_action_from_trial_start(end); 20+state_sbj.action_history(state_sbj.index); rwd_condition]]; % event save
                    end
                end
            end
            
            %% [3] moving to the next state
            [state_sbj map_sbj]=StateSpace(state_sbj,map_sbj);  % map&state index ++
            
            
            

        end % end of each choice
        
        
        %% [4] terminal state: display reward
        current_state=state_sbj.state_history(state_sbj.index);
        
        % (0) display fixation mark and display off during the jittered interval
        Screen('FillRect',wPtr,BackgroundColor_Cue_page);
        DrawFormattedText(wPtr, '.', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
        Screen(wPtr, 'Flip');
        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
            (GetSecs-session_clock_start); (GetSecs-block_clock_start); (GetSecs-trial_clock_start); 0.5; rwd_condition]]; % event save
        sec_stim_interval0=rand*(max(sec_stim_interval)-min(sec_stim_interval))+min(sec_stim_interval);
        WaitSecs(sec_stim_interval0);
        
        % (1) add state image
        Screen('FillRect',wPtr,BackgroundColor_Cue_page);
        input_stim = Screen('MakeTexture', wPtr, img_set{1,current_state});
        xpos = round(screenWidth/2);    ypos = round(screenHeight/2);
        sx=floor(IMAGE_SIZE(1)*disp_scale);       sy=floor(IMAGE_SIZE(2)*disp_scale);
        destrect=[xpos-sx/2,ypos-sy/2,xpos+sx/2,ypos+sy/2];
        Screen('DrawTexture', wPtr, input_stim,[],destrect);

        % (2) outcome message read
        outcome_state=state_sbj.state_history(3);
        original_rwd=map_sbj.reward(outcome_state);
        actual_rwd=state_sbj.reward_history(3);
        % find the corresponding outcome image
        type_outcome_final=-99;
        for oo=1:1:size(map_sbj.goal_state_index,2)
            if(length(find(map_sbj.goal_state_index{1,oo}==outcome_state))>0)
                type_outcome_final=oo;
            end
        end
        do_show_coin_img=0;
        if(original_rwd==actual_rwd) % earn case : get what the outcome state is supposed to provide.
            case_earn=1;
            msg_out=sprintf('%d added to your total.',original_rwd);            
            if(type_outcome_final~=-99) % for non-zero outcome states
                input_stim2 = Screen('MakeTexture', wPtr, img_set_reward_type{1,type_outcome_final});            
                do_show_coin_img=1;
            end
        else % devalued case - lost case
            case_earn=-1;
            msg_out=sprintf('Sorry, you do not get %d this time.',original_rwd);
            do_show_coin_img=0;
            if(type_outcome_final~=-99) % for non-zero outcome states
                input_stim2 = Screen('MakeTexture', wPtr, img_set_reward_type{1,type_outcome_final});
                do_show_coin_img=1;
            end
        end

        % (3) add outcome message
        sx2=floor(COIN_IMG_SIZE(1)*disp_scale_outcome_msg);       sy2=floor(COIN_IMG_SIZE(2)*disp_scale_outcome_msg);
        xpos2=xpos-200; ypos2=destrect_state(4)+sy2/2+10;
        destrect2=[xpos2-sx2/2,ypos2-sy2/2,xpos2+sx2/2,ypos2+sy2/2];
        if(do_show_coin_img==1) % for non-zero outcome states
            Screen('DrawTexture', wPtr, input_stim2,[],destrect2); % coin img
        end
        oldTextSize=Screen('TextSize', wPtr ,text_size_rwd_msg);
        DrawFormattedText(wPtr, msg_out, destrect2(3)+10 , ypos2-10, COLOR_FIXATION_MARK);
        oldTextSize=Screen('TextSize', wPtr ,oldTextSize);

        % (4) display on
        Screen(wPtr, 'Flip'); % display on
        if(outcome_state~=9)
            value=outcome_state+0.1*case_earn;
        else
            value=outcome_state;
        end
        onset_state_from_trial_start=[onset_state_from_trial_start (GetSecs-trial_clock_start)];
        HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
            (GetSecs-session_clock_start); (GetSecs-block_clock_start); onset_state_from_trial_start(end); value; rwd_condition]]; % event save
       
        if(strcmp(session_opt,'pre')==1) % pre-session
            WaitSecs(1.5);
        end
        if(strcmp(session_opt,'fmri')==1) % fmri-session
            WaitSecs(sec_reward_display);
        end

        %% [5] ITI - show blank screen
        if((trial<max(trial_set))&&(block<max(Tot_block)))
            Screen('FillRect',wPtr,BackgroundColor_Cue_page);
            DrawFormattedText(wPtr, 'Now ready for the next trial...', 'center', 'center', COLOR_FIXATION_MARK); % add 'o' mark at the click pt.
            Screen(wPtr, 'Flip');
            HIST_event_info0=[HIST_event_info0 [block; trial_in_block; state_sbj.index; ...
                (GetSecs-session_clock_start); (GetSecs-block_clock_start); (GetSecs-trial_clock_start); 0.5; rwd_condition]]; % event save
            sec_trial_interval0=rand*(max(sec_trial_interval)-min(sec_trial_interval))+min(sec_trial_interval);
            WaitSecs(sec_trial_interval0);
        end


        %% Update HIST_event_info (overwrite at each block)
        HIST_event_info{1,index_num}=HIST_event_info0;
        
        %% update behavior matrix !@#$%
        if(isempty(HIST_behavior_info0))
            acc_rwd=0;
        else
            acc_rwd=HIST_behavior_info0(end,17);
        end
        mat_update=[block, trial_in_block, block_condition, ...
            state_sbj.state_history(1), state_sbj.state_history(2), state_sbj.state_history(3), ...
            state_sbj.action_history(1), state_sbj.action_history(2), ...
            state_sbj.RT(1), state_sbj.RT(2), ...
            onset_state_from_trial_start(1), onset_state_from_trial_start(2), onset_state_from_trial_start(3), ...
            onset_action_from_trial_start(1), onset_action_from_trial_start(2), ...
            state_sbj.reward_history(end), acc_rwd+state_sbj.reward_history(end), ...
            rwd_condition'];
        HIST_behavior_info0=[HIST_behavior_info0; mat_update];
        
    end % end of each trial


    % behavior matrix update
    HIST_behavior_info{1,index_num}=HIST_behavior_info0;

    %% save the (updated) image usage matrix (overwriting)
    file_imgind_sv_name=[EXP_NAME '_info.mat'];
    file_name_sv=[pwd '\result_save\' file_imgind_sv_name];
    save(file_name_sv,'HIST_event_info','HIST_event_info_Tag','HIST_behavior_info','HIST_behavior_info_Tag','HIST_block_condition','HIST_block_condition_Tag');


end % end of each block

        
%% Ending message
str_end=sprintf('- Our experiments is over. Press any key to quit. -');
DrawFormattedText(wPtr, str_end, 'center', 'center');
Screen(wPtr, 'Flip');
KbWait; % temporarily disabled for test APR 21
% take a snapshot
if(DO_TAKE_SNAPSHOT==1)
    snapshot=Screen(wPtr, 'GetImage', [1, 1, floor(1.0*SCREEN_RESOLUTION(1)), floor(1.0*SCREEN_RESOLUTION(2))]);
    imageArray=[imageArray; {snapshot}];
end


%% save snapshots : CAUTION HEAVY PROCESS - might take a minute.
if(DO_TAKE_SNAPSHOT==1)
    for j=1:1:size(imageArray,1)
        str=sprintf('snapshot_dispay_exp_%03d.png',j);
        imwrite(imageArray{j},['snapshot\' str],'png');
    end
end


% behavior matrix update
HIST_behavior_info{1,index_num}=HIST_behavior_info0;


%% save the (updated) image usage matrix (overwriting)
file_imgind_sv_name=[EXP_NAME '_info.mat'];
file_name_sv=[pwd '\result_save\' file_imgind_sv_name];
save(file_name_sv,'HIST_event_info','HIST_event_info_Tag','HIST_behavior_info','HIST_behavior_info_Tag','HIST_block_condition','HIST_block_condition_Tag');
if(strcmp(session_opt,'pre')==1)
    if(index_num==1) % for the first session per each subject
        save(file_name_sv,'HIST_event_info','HIST_event_info_Tag','HIST_behavior_info','HIST_behavior_info_Tag','HIST_block_condition','HIST_block_condition_Tag','Info_seed');
    end
end

%% save all variables
file_sv_name=[EXP_NAME sprintf('_%d.mat',index_num)];
file_name_sv=[pwd '\result_save\' file_sv_name];
save(file_name_sv,'*');


%% session end sound
sec_dur_sound=1;
Beeper('med', 0.4, sec_dur_sound); WaitSecs(sec_dur_sound)
Beeper('high', 0.4, sec_dur_sound)



%% finish all
Screen('CloseAll');
clear mex
% clear Screen

disp('########################################################')
str_end1=sprintf('### session%d is done ############################',index_num);
disp(str_end1);
str_end2=sprintf('### next session = %d ############################',index_num+1);
disp(str_end2);
disp('########################################################')

% display the number of response failure
missed_count = length(find(HIST_event_info{1,index_num}(7,:)==-99));
disp(sprintf('- # of response failure in this session = %d. (will be penalized) ',missed_count));

output_info=1;

end


