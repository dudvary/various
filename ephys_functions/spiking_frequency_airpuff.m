%% Analyze Spiking of Air Puff recordings
% Author: Daniel Udvary
% Version Date: 1st July 2020
% Input:
%   - filepath
%   - baseline_time = 60000; % ms (optional)
%   - airpuff_duration = 700; % ms (optional)
%   - ISI = 2.5*1000; % ms (2.5s) (optional)
%   - trials_num = 30; (optional)
%   - binsz = 5; % ms (optional)
%   - r = [-200 1100]; % range given airpuff (t=0) (optional)
% Output:
%   - Plots PSTH over 30 trials in range [-200ms to 1100ms] (airpuff at t=0)
%   - .csv file containing total number of spikes over all trials in given time range (stored in same folder)
%   - Displays Baseline, Airpuff, and AfterAirpuff Frequency
% Bug fixes:
%   - 1st July 2020: Multiple spikes were missed if in same bin
function spiking_frequency_airpuff(filepath,baseline_time,airpuff_duration,ISI,trials_num,binsz,r)
        
    % [X TrialNum X Spiketimes [ms]] 
    data = ReadCluster(filepath);
    [~,filename,~] = fileparts(filepath);   
    
    switch nargin
        case 1
            baseline_time = 60000; % ms
            airpuff_duration = 700; % ms
            ISI = 2.5*1000; % ms (2.5s)
            trials_num = 30; 
            binsz = 5; % ms
            r = [-200 1100]; % range given airpuff (t=0)
        case 2
            airpuff_duration = 700; % ms
            ISI = 2.5*1000; % ms (2.5s)
            trials_num = 30; 
            binsz = 5; % ms
            r = [-200 1100]; % range given airpuff (t=0)
        case 3
            ISI = 2.5*1000; % ms (2.5s)
            trials_num = 30; 
            binsz = 5; % ms
            r = [-200 1100]; % range given airpuff (t=0)
        case 4
            trials_num = 30; 
            binsz = 5; % ms
            r = [-200 1100]; % range given airpuff (t=0)
        case 5
            binsz = 5; % ms
            r = [-200 1100]; % range given airpuff (t=0)
        case 6
            r = [-200 1100]; % range given airpuff (t=0)
    end
   
    % Check data
    if (length(unique(data(:,2)))>1)
        warning(['More than one trial detected! Number of trials detected: ' ...
            num2str(length(unique(data(:,2))))]); 
    end
    
    mxSpikeTime = max(data(:,4)); 
    if (mxSpikeTime < baseline_time)
        warning(['Last Spikes only detected during Baseline Activity' ...
            'Last Spike at ' num2str(mxSpikeTime) 'ms']);
    elseif (ceil((mxSpikeTime-baseline_time)/(airpuff_duration+ISI))+1 < trials_num)
        tmp = ceil((mxSpikeTime-baseline_time)/(airpuff_duration+ISI))+1;
        warning(['Last Spike only detected in trial ' num2str(tmp) ' of ' ...
            num2str(trials_num) '. Last Spike at ' num2str(mxSpikeTime) 'ms']);
    end
        
    fprintf('Settings of Protocol:\n');
    fprintf('  Duration of Airpuff: %dms\n',airpuff_duration);
    fprintf('  Duration of ISI: %dms\n',ISI);
    fprintf('  Number of trials: %d\n',trials_num);
    fprintf('\n');

    baseline_freq = sum(data(2:end,4)<baseline_time)/(baseline_time/1000); 
    fprintf('Baseline Frequency: %fHz\n',baseline_freq);
  
    binSpikesTrial = zeros(trials_num,sum(abs(r))/binsz);
    for t = 1:trials_num
        start_t = baseline_time+(t-1)*(ISI+airpuff_duration) + r(1); 
        end_t = start_t + (r(2)-r(1)); 
        tmp = (data(:,4)<end_t)+(data(:,4)>=start_t);
        tmp_times = data(tmp==2,4)-start_t; 
        tmp_binID = floor(tmp_times./binsz)+1;
        
        for tt = tmp_binID'
            binSpikesTrial(t,tt) = binSpikesTrial(t,tt)+1; 
        end
    end
    binSpikes = sum(binSpikesTrial);
    
    % Baseline Frequency
    idxBase = floor(abs(r(1))/binsz); 
    freqBaseline = sum(binSpikes(1:idxBase)) / trials_num * 1000/abs(r(1)); 
    fprintf('  Spontenous Frequency [%d 0]: %fHz\n',r(1),freqBaseline);
    
    % Airpuff Frequency
    idxAirpuff = floor(airpuff_duration/binsz)+idxBase; 
    freqAirpuff = sum(binSpikes(idxBase+1:idxAirpuff)) / trials_num * 1000/airpuff_duration; 
    fprintf('  Airpuff Frequency [0 %d]: %fHz\n',airpuff_duration,freqAirpuff);
    
    % After-Airpuff Frequency
    idxAfterAirpuff = floor((r(2)-airpuff_duration)/binsz)+idxAirpuff; 
    freqAfterAirpuff = sum(binSpikes(idxAirpuff+1:idxAfterAirpuff)) / trials_num * 1000/(r(2)-airpuff_duration); 
    fprintf('  After-Airpuff Frequency [%d %d]: %fHz\n',airpuff_duration,r(2),freqAfterAirpuff);
        
    %% Write CSV file
    filenameCSV = filepath(1:end-9); 
    if ~(strcmp(filenameCSV(end-3:end), '.csv'))
        filenameCSV = [filenameCSV '.csv'];
    end
             
    fid = fopen(filenameCSV, 'w');
    x = [0:binsz:sum(abs(r))];
    x(end) = []; 
    x = x+binsz/2; 
    
    fprintf(fid, '%s, %d trials, %d to %dms baseline, %d to %dms airpuff, %d to %dms afterairpuff \n', ...
        filename,trials_num, 0, abs(r(1)), abs(r(1)), abs(r(1))+airpuff_duration, ...
        abs(r(1))+airpuff_duration, sum(abs(r)));
    fprintf(fid, 'binEdge1[ms],binEdge2[ms],number_of_spikes\n');
    
    for i = 1:length(binSpikes)
        fprintf(fid, '%f,%f,%f\n', x(i)-binsz/2, x(i)+binsz/2, binSpikes(i));
    end
    fclose(fid);
    
    %% Plot PSTH
    binSpikesNorm = binSpikes./trials_num; 
    f1 = figure();
    whitebg('white');
    x = [r(1):binsz:r(2)];
    x(end) = []; 
    x = x+binsz/2; 
    bar(x,binSpikesNorm,'k'); 
    hold on;
    plot([0 0],[0 max(binSpikesNorm)*1.1],'r--');
    hold on;
    plot([airpuff_duration airpuff_duration],[0 max(binSpikesNorm)*1.1],'r--');
    ylabel(['APs per ' num2str(binsz) 'ms']);
    xlabel('time [ms]');
    set(gca,'XLim',[r(1) r(2)],'YMinorTick','on','TickDir','out', ...
     'Box','off','XMinorTick','on','YLim',[0 max(max(binSpikesNorm)*1.1,0.01)], ...
     'XTick',[r(1):100:r(2)])
 
    set(f1,'PaperPositionMode', 'auto','Position',[0 0 1500 1000]);
    print(f1,'-depsc2','-r600', [filepath(1:end-9) '_PSTH']); 
    print(f1,'-dpng','-r600',[filepath(1:end-9) '_PSTH']);
 
    %% Plot Raster
    f2 = figure();
    for t = 1:trials_num
         timetmp = (find(binSpikesTrial(t,:)>0)-1)*binsz+binsz/2 + r(1); 
         for tt = timetmp
            plot([tt tt],[t-0.3 t+0.3],'k-','LineWidth',2);
            hold on;
         end
         hold on; 
    end
    hold on;
    plot([0 0],[0 trials_num+1],'r--');
    hold on;
    plot([airpuff_duration airpuff_duration],[0 trials_num+1],'r--');
    
    xlabel('time [ms]');
    ylabel('trial #');
    set(gca,'XLim',[r(1) r(2)],'YMinorTick','on','TickDir','out', ...
        'Box','off','XMinorTick','on','YLim',[0.5 trials_num+0.5], ...
        'XTick',[r(1):100:r(2)])
    
    set(f2,'PaperPositionMode', 'auto','Position',[0 0 1500 1000]);
    print(f2,'-depsc2','-r600', [filepath(1:end-9) '_Raster']); 
    print(f2,'-dpng','-r600',[filepath(1:end-9) '_Raster']);
end