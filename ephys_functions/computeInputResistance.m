function RI = computeInputResistance(voltageTrace,currentTrace)
% Compute Somatic InputResistance: Slope of fitted line to Voltage-Current
% Curve
% Input:
% - voltageTrace: [NumTimeSteps x CurrentPulses]
%           Voltage Traces for each Current Pulse (in milliVolt)
% - currentTrace: [NumTimeSteps x CurrentPulses]
%           Current Trace for each Current Pulse (in pikoAmpere)
% Output:
% - RI: Somatic Input Resistance (in MegaOhm)

% DEFINITION: (L6INs)
% somatic input resistance was defined as the slope of a current / voltage 
% curve based on 20 decreasing rectangular current pulses (starting near 
% the previously defined AP firing threshold and typically reaching 
% baseline or slightly hyperpolarized membrane potential with the last
% pulse).

    if size(currentTrace,1) ~= size(voltageTrace,1) || size(currentTrace,2) ~= size(voltageTrace,2)
       error('Number of Voltage Traces does not match number of Current Traces');  
    end
    
    if size(currentTrace,2) ~= 20
       warning('Not 20 Current Pulses applied (as done for L6INs)'); 
    end

    % Find indices when CurrentPulse was applied
    % current trace should be ~= 0
    tspanCurrentPulse = [];
    k = 1;
    while isempty(tspanCurrentPulse)
        tspanCurrentPulse = find(currentTrace(:,k)~=0);
        k = k+1;
    end 

    if length(tspanCurrentPulse)~=2500
        warning(['No proper Current Pulse found']);
    end

    % Compute median Voltage Potential
    medVoltageRi = median(voltageTrace(tspanCurrentPulse,:)); 

    % Current is in pikoOhm (E-12)
    xData = currentTrace(tspanCurrentPulse(1),:).*1e-12;
    % Voltage is in miliVolt (E-3)
    yData = medVoltageRi.*1e-3;

    % Do polyfit1 a*x+b
    warning('off');    
    p = polyfit(xData,yData,1);
    warning('on');
    
    % Slope is InputResistance R = U/I
    RI = p(1)/1e6; % in MOhm    
end