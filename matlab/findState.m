function [state, port] = findState(state_prefix, states)
    state = states(strncmpi(state_prefix, states, numel(state_prefix)));
    if isempty(state)
        state = state_prefix;
        port = '';
    else
        state = state{1};
    % This only works if your state has one poke name at the end.
    % e.g. reward_state_BotC
        port = state(numel(state_prefix)+2:end);
    end