function [state_times] = getStateChange(pevents, state_name, varargin)
% [e_times] = getStateChange(pevents, StateName, getin)
%
% Inputs:
% pevents          A array of parsed events
%
% Optional Input
% change_types    ['first-in'] which state change events to include
%                   'first-in' returns just the first entry of that state per trial
%                   'first-out' returns just the first exit of that state per trial
%                   'last-in' returns just the last entry of that state per trial
%                   'last-out' returns just the last exit of that state per trial

% prefix          [false] if true then your state_name will match against
%                   states names that have an extra _suffix.  e.g. state1
%                   will match state1_BotC
%
%
% Output:
% event_times    a list of the times relative to the first trial in the input

if nargin==0
    help(mfilename)
    return
end

inpd = @utils.inputordefault;

change_types = inpd('change_types','first-in',varargin);
prefix = inpd('prefix',false, varargin);

session_start = pevents(1).StartTime;
state_times = nan(numel(pevents),1);

if prefix
    switch change_types
        case 'first-in'
            sx_fcn = @(tx)int_fcn_first(pevents(tx),state_name);
        case 'first-out'
            sx_fcn = @(tx)int_fcn_last(pevents(tx),state_name);
        case 'last-in'
            sx_fcn = @(tx)int_lcn_first(pevents(tx),state_name);
        case 'last-out'
            sx_fcn = @(tx)int_lcn_last(pevents(tx),state_name);
    end
else
    switch change_types
        case 'first-in'
            sx_fcn = @(tx)pevents(tx).States.(state_name)(1,1);
        case 'first-out'
            sx_fcn = @(tx)pevents(tx).States.(state_name)(1,2);
        case 'last-in'
            sx_fcn = @(tx)pevents(tx).States.(state_name)(end,1);
        case 'last-out'
            sx_fcn = @(tx)pevents(tx).States.(state_name)(end,2);
    end
end
    % this uses prefix matching.
for tx = 1:numel(pevents)
    try
        [t_events] = sx_fcn(tx);
        state_times(tx) = t_events + pevents(tx).StartTime;
    catch me%  should only catch MATLAB:nonExistentField
        if ~strcmp(me.identifier,'MATLAB:nonExistentField')
            rethrow(me)
        end
    end

end
            
state_times = state_times - session_start;


end

function et = int_fcn_first(pes, sn)
    full_state_name = findState(sn, fields(pes.States));
    et = pes.States.(full_state_name)(1,1);

end

function et = int_fcn_last(pes, sn)
    full_state_name = findState(sn, fields(pes.States));
    et = pes.States.(full_state_name)(1,2);
end
function et = int_lcn_first(pes, sn)
    full_state_name = findState(sn, fields(pes.States));
    et = pes.States.(full_state_name)(end,1);

end
function et = int_lcn_last(pes, sn)
    full_state_name = findState(sn, fields(pes.States));
    et = pes.States.(full_state_name)(end,2);

end
