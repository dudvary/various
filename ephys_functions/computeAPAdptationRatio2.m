function [APAdaptationRatio,ISI] = computeAPAdptationRatio2(tspan,spikeIndices)
% computeAPAdaptationRatio2
% Input:
% - tspan: [1 x NumTimeSteps]
% - spikeIndices: [1 x NumSpikes] indices of spikes in tspan
% Output:
% - APAdaptationRatio: ratio of last and first point of a monoexponential 
% fit to the spiking frequency.
% - ISI: Inter-Spike-Interval

% DEFINITION (L4NFS):
% AP adaptation ratio was defined as the ratio b/w the avg of 1st three
% and last three spikes in a 10-spike train. 

    if min(spikeIndices)<1 || max(spikeIndices)>length(tspan)
       error('Indices in spikeIndices exceeds Number of Time Steps'); 
    end

    % Inter-Spike-Interval
    ISI = diff(tspan(spikeIndices));    
    % APAdaptationRatio
    APAdaptationRatio = mean(ISI(1:3))./mean(ISI(end-2:end));
end

