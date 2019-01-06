%% pre-training 
Screen('Preference', 'SkipSyncTests', 1);
% 1. insert a subject's name.
subcode = 'kdj';
% 2. pretraining 
SIMUL_cogload_behavior(subcode, 1, 'pre');

%% main
% session 1
SIMUL_cogload_behavior(subcode, 1, 'fmri');
% session 2
SIMUL_cogload_behavior(subcode, 2, 'fmri');
% session 3
SIMUL_cogload_behavior(subcode, 3, 'fmri');
% session 4
SIMUL_cogload_behavior(subcode, 4, 'fmri');

% KEY_L1=49; %1
% KEY_L2=50; %2
% KEY_R2=51; %3
% KEY_R1=52; %4
% KEY_Y=89; %'y'
% KEY_N=78; %'n'
% KEY_Q=81; %'q'
% KEY_T=84; % 't'
