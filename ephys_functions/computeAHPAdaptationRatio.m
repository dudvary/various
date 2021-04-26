function [AHPAdaptationRatio,AHP] = computeAHPAdaptationRatio(voltageTrace,spikeIndices)
%computeAHPAdaptationRatio 
% Input:
% - voltageTrace: [NumTimeSteps x 1]
% - spikeIndicies: [1 x NumSpikes] indices of spikes in tspan
%       (optional)
% Output:
% - AHPAdaptationRatio: see Definition
% - AHP: AfterHyperpolarizationPotential (minimum between to spikes) [1 x
%           NumSpikes-1]

% DEFINITION (L6INs):
% The minimum between two spikes was plotted as after hyperpolarization (AHP), 
% The ratio of the last and the first point of a monoexponential fit to 
% that plot was defined as AHP adaptation ratio 

    if nargin<2
        % Find spikes
        Vmod = voltageTrace;
        Vmod(Vmod<0) = 0; 
        [~,spikeIndices] = findpeaks(Vmod); 
    end
    if nargin==2 && (min(spikeIndices)<1 || max(spikeIndices)>length(voltageTrace))
       error('Indices in spikeIndices exceeds Number of Time Steps'); 
    end
    
    % Monoexponetial Fit
    ft = fittype('a*exp(b*x)+c'); 

    % AfterHyperpolarizationPotential 
    % Compute Minimum between two spikes
    AHP = nan(1,length(spikeIndices)-1);
    for i = 1:length(spikeIndices)-1
        AHP(i) = min(voltageTrace(spikeIndices(i):spikeIndices(i+1)));
    end

    % Fit Data to monoexponential 
    xData = 1:length(AHP);
    try
        f = fit(xData',AHP',ft,'StartPoint',[0 0 0]);
    catch
        warning('Fitting failed in computeAHPAdaptationRatio'); 
        AHPAdaptationRatio = nan;
        return;
    end
    c = coeffvalues(f);
    yFit = c(1).*exp(c(2).*xData)+c(3); 

    AHPAdaptationRatio = yFit(end)/yFit(1);
    AHPAdaptationRatio = AHPAdaptationRatio*sign(c(2)); 
end

