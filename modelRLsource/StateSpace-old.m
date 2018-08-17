function [myState]=StateSpace(myState,myMap)
% input: myState (with current ".state_history" and ".action_history"
% parameter: myMap.action1, myMap.action2, myMap.reward
% output: myState (with updated ".index", ".state_history" and ".reward_history"

%% readout
current_state=myState.state_history(myState.index);
current_action=myState.action_history(myState.index);


if(myMap.IsTerminal(current_state)~=1)
    
    %% clock+1    
    myState.index=myState.index+1;
    
    %% state transition
    prob_mat=myMap.action(1,current_action).prob(current_state,:); %(ex) [... 0 0.3 0 0.7]
    state_mat=myMap.action(1,current_action).connection(current_state,:); % (ex) [... 0 1 0 1]
    [tmp state_cand]=find(state_mat==1); % a set of candidate next-states
    if(rand<=prob_mat(state_cand(1)))
        % state transition
        myState.state_history(myState.index)=state_cand(1);
        % reward
        myState.reward_history(myState.index)=myMap.reward(state_cand(1));
    else
        % state transition
        myState.state_history(myState.index)=state_cand(2);
        % reward
        myState.reward_history(myState.index)=myMap.reward(state_cand(2));
    end
    
    %% filling in SARSA matrix - (r,s')
    myState.SARSA(3:4)=[myState.reward_history(myState.index) myState.state_history(myState.index)];
    
    %% If the new state is terminal,
    if(myMap.IsTerminal(myState.state_history(myState.index))==1)
        myState.JobComplete=1;
    end
        
end


end