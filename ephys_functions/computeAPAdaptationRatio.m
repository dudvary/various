function APAdaptationRatio = computeAPAdaptationRatio(tspan,spikeIndices)
% computeAPAdaptationRatio
% Input:
% - tspan: [1 x NumTimeSteps]
% - spikeIndices: [1 x NumSpikes] indices of spikes in tspan
% Output:
% - APAdaptationRatio: ratio of last and first point of a monoexponential 
% fit to the spiking frequency. 

% DEFINITION (L6IN):
% AP adaptation ratio was defined as the ratio of the last and first 
% point of a monoexponential fit to the spiking frequency. 

    if min(spikeIndices)<1 || max(spikeIndices)>length(tspan)
       error('Indices in spikeIndices exceeds Number of Time Steps'); 
    end

    % Monoexponetial Fit
    ft = fittype('a*exp(b*x)+c'); 

    % Spike Frequency
    ISI = diff(tspan(spikeIndices));    % time between spikes
    spikeFrequency = 1./ISI;            % SpikeFrequency
    
    % Fit Data to monoexponential 
    xData = 1:length(spikeFrequency);
    
    try     
        f = fit(xData',spikeFrequency',ft,'StartPoint',[0 0 0]);
    catch
        warning('Fitting failed in computeAPAdaptationRatio'); 
        APAdaptationRatio = nan;
        return;
    end
    c = coeffvalues(f);
    yFit = c(1).*exp(c(2).*xData)+c(3); 

    APAdaptationRatio = yFit(end)/yFit(1);
    APAdaptationRatio = APAdaptationRatio*sign(c(2));
    
end