%% [fitResults] = getMLEfits(YCount,YSample,XCount,XSample)
% Does MLE fit for discrete distributions:
%   'Binomial','Poisson','Negative Binomial'
% Does MLE fit for continous distributions;
%   'Half Normal','Normal','Gamma','Exponential'
% Input: 
% - YCount: Count data (i.e., connection probabilities mapped on integers
%       from 0 to 100)
% - YSample: Connection probabilities
% - XCount: Bins for count data (0 to 100)
% - XSample: Bins for connection probabilities (0 to 1)
% Output:
% - fitResults: Structure containing fit results
%     fitResults.yFit: fitted y values
%     fitResults.y: original values (empty by default)
%     fitResults.x: bins
%     fitResults.pd: probability distribution class (by Matlab)
%           contains a lot more information
%     fitResults.distName: Name of distribution 
function [fitResults] = getMLEfits(YCount,YSample,XCount,XSample)

    ListOfDiscreteDist = {'Binomial','Poisson','Negative Binomial'};
    ListOfContinousDist = {'Half Normal','Normal','Gamma','Exponential'};
    
    N = length(ListOfDiscreteDist)+length(ListOfContinousDist); 
        
    fitResults.yFit = cell(1,N);
    fitResults.y = cell(1,N);
    fitResults.x = cell(1,N);
    fitResults.pd = cell(1,N);
    fitResults.distName = cell(1,N); 
    
    c = 1; 
    for i = 1:length(ListOfDiscreteDist)        
        
        if strcmp(ListOfDiscreteDist{i},'Binomial')
            nTrials = 100; 
            pd = fitdist(YSample.*nTrials,ListOfDiscreteDist{i}, ...
                    'ntrials',nTrials);
        else
            pd = fitdist(YCount,ListOfDiscreteDist{i});
        end
        
        fitResults.pd{c} = pd; 
        fitResults.x{c} = XCount; 
        fitResults.y{c} = [];
        fitResults.yFit{c} = pdf(pd,XCount); % all ready as probabilities
        fitResults.distName{c} = ListOfDiscreteDist{i};
        c = c+1; 
    end
    
    for i = 1:length(ListOfContinousDist)
        
        ytmp = YSample;
        xtmp = XSample;

        if strcmp(ListOfContinousDist{i},'Gamma')
            ytmp(ytmp==0) = eps;
            xtmp(xtmp==0) = xtmp(2)/2; 
        end
        
        pd = fitdist(ytmp,ListOfContinousDist{i});
        yFit = pdf(pd,xtmp); 
        % yFit does not represent probabilities
        % therefore normalize by sum and then subtract area outside of 0
        % and 1
        b = [cdf(pd,0) 1-cdf(pd,1)]; 
        yFit = yFit./sum(yFit);
        yFit = yFit.*(1-sum(b));   
        
        fitResults.pd{c} = pd; 
        fitResults.x{c} = xtmp; 
        fitResults.y{c} = [];
        fitResults.yFit{c} = yFit; % all ready as probabilities
        fitResults.distName{c} = ListOfContinousDist{i};

        c = c+1; 
    end
end