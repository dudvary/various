%% Compute probability of Chain Motifs across different motif sizes
% and compare probability with random/avg network model (only constraint by
% average connection probability)
% Chain motif: every node but two has one incoming edge and one outcoming edge. 
%   one node has only one incoming edge (sink), 
%   and one other only one outcoming edge (source).
clear all
close all
clc

% Please modify
filename = ''; % path/to/matrix containing connection probabiltiy values
% should be square!
outputPath = ''; % path where figures should be saved!
load(filename,'pMatrix'); 

numTrials = 10000;
numChainMotifs = 1000; 
numNodesList = [3:1:10]; 
n = size(pMatrix,1); 

%%
rng(154581); 
for numNodes = numNodesList
    
    edgeMatrix = generateRandomChainMotif(numNodes,numChainMotifs);
    pSelect = nan(numTrials,numNodes,numNodes); 
    idDiagonal = 1:numNodes+1:(numNodes^2); % id of diagonal values

    for j = 1:numTrials
        idx = randperm(n,numNodes);
        p = pMatrix(idx,idx);
        p(idDiagonal) = nan; 
        pSelect(j,:,:) = p; 
    end

    p_avg = nanmean(pSelect(:));
    p_sd = nanstd(pSelect(:)); 
    
    [edgeStats] = extractEdgeMotifProbabilityGivenChainMotif( ...
                                edgeMatrix,pSelect,p_avg);     
    edgeStats.comment = 'chain'; 
    edgeStats.p_sd = p_sd; 

    % SAVE DATA
    save([outputPath 'ChainMotif_NumTrials_' num2str(numTrials) ...
        '_NumCells_' num2str(numNodes) '.mat'],'edgeStats');
    fprintf('SAVED for %d\n',numNodes); 
end