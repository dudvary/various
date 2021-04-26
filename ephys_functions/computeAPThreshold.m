function APthreshold = computeAPThreshold(voltageTrace,tspan,tstim)
% compute Action Potential Threshold
% Input:
% - voltageTrace [NumTimeSteps x NumVoltageTraces] 
%       has to be sorted according to current pulse amplitude
% - tspan [1 x NumTimeSteps] [in msec]
% - tsim [1 x 2] [in msec] 
%       start and ending time of current pulse (usually 100 and 600)
% Output:
% - APthreshold
%       Average of median Voltage of last trace with no spike and first
%       trace with a spike elicted by increasing rectengular current pulses

% tested for L6 INs

% DEFINITION: (L6INs)
% "The AP threshold was calculated as the average of the median membrane 
% potential during the current pulse of the last trace without and the 
% first trace with an AP."

    if length(tstim)~=2
       error('tstim has to be 2-array vector');
    end

    if size(voltageTrace,2) == length(tspan)
        voltageTrace = voltageTrace';
    end

    if size(voltageTrace,1) ~= length(tspan)
       error('Length of Vthreshold does not match length of tspan');
    end
    
    if size(voltageTrace,2) < 2
       error('Vthreshold contains only one Voltage Trace, two needed!');
    end
    
    if tstim(1)<tspan(1) || tstim(2) > tspan(end)
        error('Time of Stimulus does not match TimeTrace');
    end

    APthreshold = nan; 

    % Find trace with the first spike
    tmpDiff = diff(max(voltageTrace));
    [~,idx] = max(tmpDiff);
    idx = idx+1;

    if isempty(idx)
       warning(['No Spike in Vthreshold detected!']);
       return;
    end

    medVoltage1 = median(voltageTrace(tspan>=tstim(1) & tspan<=tstim(2),idx-1)); 
    medVoltage2 = median(voltageTrace(tspan>=tstim(1) & tspan<=tstim(2),idx)); 
    APthreshold = mean([medVoltage1 medVoltage2]);
end