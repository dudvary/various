%% Compute probability of fully recurrent network motif across different motif sizes
% and compare probability with random/avg network model (only constraint by
% average connection probability)
clear all
close all
clc

% Please modify
filename = ''; % path/to/matrix containing connection probabiltiy values
% should be square!
outputPath = ''; % path where figures should be saved!

%%
load(filename,'pMatrix'); 
numTrials = 1e7;
numNodesListFullyRec = [3:1:10]; 
n = size(pMatrix,1); 

%%
rng(154581); 
pRec = nan(size(numNodesListFullyRec)); 
pRecRandom = nan(size(numNodesListFullyRec)); 
devFullyRec = nan(size(numNodesListFullyRec));
p_avg = nan(size(numNodesListFullyRec)); 
p_sd = nan(size(numNodesListFullyRec));

for numNodes = numNodesListFullyRec

    idDiagonal = 1:numNodes+1:(numNodes^2); % id of diagonal values
    idNonDiagonal = setdiff(1:numNodes^2,idDiagonal); 
    pRec_tmp = nan(1,numTrials); 
    
    pSelect = nan(numTrials,numNodes,numNodes); 
    
    for j = 1:numTrials
        idx = randperm(n,numNodes);
        p = pMatrix(idx,idx);
        pRec_tmp(j) = prod(p(idNonDiagonal));  
        p(idDiagonal) = nan;
        pSelect(j,:,:) = p;
    end
    
    p_avg(numNodes==numNodesListFullyRec) = nanmean(pSelect(:));
    p_sd(numNodes==numNodesListFullyRec) = nanstd(pSelect(:)); 
    
    pRec(numNodes==numNodesListFullyRec) = mean(pRec_tmp); 
    pRecRandom(numNodes==numNodesListFullyRec) = ...
            p_avg(numNodes==numNodesListFullyRec)^numel(idNonDiagonal);     

    devFullyRec(numNodes==numNodesListFullyRec) = ...
                    pRec(numNodes==numNodesListFullyRec)./ ...
                    pRecRandom(numNodes==numNodesListFullyRec);
    
    fprintf('%d,%.2e,%.2e,%.2e\n',numNodes, ...
            pRec(numNodes==numNodesListFullyRec), ...
            pRecRandom(numNodes==numNodesListFullyRec), ...
            devFullyRec(numNodes==numNodesListFullyRec));
end

%%
figure(1);
clf; 
plot(numNodesListFullyRec,devFullyRec,'k.-'); 
hold on;
plot(numNodesListFullyRec,pRec,'r.-'); 
plot(numNodesListFullyRec,pRecRandom,'b.-'); 
plot(numNodesListFullyRec,ones(size(numNodesListFullyRec)),'k:');
set(gca,'YScale','log','Box','off','TickDir','out'); 
ylabel('Deviation');
xlabel('#Nodes'); 

% SAVE DATA
save([outputPath 'FullyRecurrentMotif_NumTrials_' num2str(numTrials) '.mat'], ...
    'pRec','pRecRandom','numNodesListFullyRec','p_avg','p_sd');
