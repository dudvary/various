function [param] = ... 
    computeEPhysParams(voltageTrace,tspan,tstim)
% Compute EPhysParameters of given voltageTrace
% Input: 
% - voltageTrace: [1 x NumTimeSteps] (in msec)
% - tspan: [1 x NumTimeSteps] (in msec)
% - tstim: [1 x 2] Start and ending time of Stimulus (in msec)
% Output:
% - param struct containing
% - APAdaptationRatio: ratio of last and first point of a monoexponential 
%       fit to the spiking frequency. 
% - AHPAdaptationRatio: see Definition
% - AHP: AfterHyperpolarizationPotential
% - APHalfWidth: HalfWidth of Action Potential

    if tstim(1)<tspan(1) || tstim(2) > tspan(end)
        error('Time of Stimulus does not match TimeTrace');
    end
    
    if length(voltageTrace) ~= length(tspan)
       error('Length of Vthreshold does not match length of tspan');
    end
    
    % Find spikes
    Vmod = voltageTrace;
    Vmod(Vmod<0) = 0; 
    [~,spikeIndices] = findpeaks(Vmod); 

    if isempty(spikeIndices)
        param.APAdaptationRatio = nan;
        param.AHP = nan; 
        param.AHPAdaptationRatio = nan;
        param.APHalfWidth = nan; 
        warning('No spikes detected in Voltage Trace');
        return; 
    elseif length(spikeIndices)<4
        param.APAdaptationRatio = nan;
        param.AHP = nan; 
        param.AHPAdaptationRatio = nan;
        param.APHalfWidth = nan; 
        warning(['Not enough Spikes in VoltageTrace detected (#spikes: ' ...
            num2str(length(spikeIndices)) ')']);
        return; 
    end
    
    % AP Adaptation Ratio
    param.APAdaptationRatio = computeAPAdaptationRatio(tspan,spikeIndices);

    % AHP adaptation ratio   
    [AHPAdaptationRatio,AHP] = computeAHPAdaptationRatio(voltageTrace,spikeIndices);
    param.AHPAdaptationRatio = AHPAdaptationRatio;
    param.AHP = AHP;

    % AP Half Width
    param.APHalfWidth = computeAPHalfWidth(voltageTrace,tspan,tstim,spikeIndices);
    
end

