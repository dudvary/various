function APHalfWidth = computeAPHalfWidth(voltageTrace,tspan,tstim,spikeIndices)
%COMPUTEAPHALFWIDTH Summary of this function goes here
%   Detailed explanation goes here
% Input: 
% - voltageTrace: [1 x NumTimeSteps] (in msec)
% - tspan: [1 x NumTimeSteps] (in msec)
% - tstim: [1 x 2] Start and ending time of Stimulus (in msec)
% - spikeIndices: [1 x NumSpikes] indices of spikes in tspan 
%       (optional)
% Output:
% - APHalfWidth: HalfWidth of Action Potential

% DEFINITION (L6INs)
% AP half-width was defined as the width of the AP at half maximal amplitude, 
% half maximal amplitude being the voltage halfway between the spike peak and 
% median voltage, averaged for all spikes of the 100 Hz trace 
    
    if nargin<4
        % Find spikes
        Vmod = voltageTrace;
        Vmod(Vmod<0) = 0; 
        [~,spikeIndices] = findpeaks(Vmod); 
    end
    if nargin==4 && (min(spikeIndices)<1 || max(spikeIndices)>length(voltageTrace))
       error('Indices in spikeIndices exceeds Number of Time Steps'); 
    end
    
    if length(tspan)~=length(voltageTrace)
        error('Length of TimeTrace does not match length of Voltage Trace');
    end
    
    if tstim(1)<tspan(1) || tstim(2) > tspan(end)
        error('Time of Stimulus does not match TimeTrace');
    end


    % Median Voltage
    medVoltage = median(voltageTrace(tspan>=tstim(1) & tspan<=tstim(2))); 
    delta_t = diff(tspan(1:2)); 

    % HalfMaximumAmplitude
    HMA = (voltageTrace(spikeIndices)-medVoltage)./2 + medVoltage;
    HW = nan(1,length(spikeIndices));

    for i = 1:length(spikeIndices)

        HMA_diff = voltageTrace - HMA(i); 

        % Find indices of values above HMA
        idx = spikeIndices(i); 

        while sum(HMA_diff(idx)<0)==0
            idx = [idx idx(end)+1];
        end
        idx(end) = [];

        while sum(HMA_diff(idx)<0)==0
            idx = [idx(1)-1 idx];
        end 
        idx(1) = []; 

        HW(i) = length(idx)*delta_t;

    end

    APHalfWidth = mean(HW); 

end

