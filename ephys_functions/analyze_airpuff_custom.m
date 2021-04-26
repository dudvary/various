%% Analyze Spiking of Air Puff recordings
% Author: Daniel Udvary
% Date: August 2015
% Plots PSTH over 30 trials in range [-200ms to 1100ms] (airpuff at t=0)
% Outputs .csv file containing total number of spikes over all trials in
% given time range (stored in same folder)
% Displays Baseline, Airpuff, and AfterAirpuff Frequency
%
% NOTE: Cluster files do not contain any trigger, thus stimulus trigger 
% is determined indirectly via baseline_time, airpuff_duration, and ISI
function analyze_airpuff_custom()

    [filename, pathname] = uigetfile('*.cluster?', 'Select a cluster file to read');
    filepath = fullfile(pathname, filename);
    baseline_time = 60000; % ms (1 min) [CAUTION: TRIGGER DEPENDENT PARAMETER]
    airpuff_duration = 700; % ms [CAUTION: TRIGGER DEPENDENT PARAMETER]
    ISI = 2.5*1000; % ms (2.5s) [CAUTION: TRIGGER DEPENDENT PARAMETER]
    trials_num = 30; 
    binsz = 0.001; 
    r = [-200 1100]; % range given airpuff (t=0) 

    spiking_frequency_airpuff(filepath,baseline_time,airpuff_duration,ISI,trials_num,binsz,r);  
end