%% Function to compute Doublet Motif Spectra of connection probability matrix
% p 
% Input: 
% - p: Connection probability matrix (N x N)
% Output:
% - result: struct
%     doubletMotifProbability: probability of unidirectional, bidirectional
%       and unconnected motif (1 x 3)
%     doubletMotifLabel: label of motif (1 x 3)
%     deviation: deviation between random/avg network and double motifs (1 x 3)
%     doubletMotifProbability_avg: probability of motifs in random/avg
%       network (1 x 3)
%     p_avg_sample: average sample probability 
%     doubletMotifProbability_avgSample: probability of motifs in random/avg
%       network (1 x 3) using same sample (should be same result)
%
function [result] = getDoubletMotifSpectra(p)

    if size(p,1)~=size(p,2)
       error('Needs to be square matrix'); 
    end
    
    if max(p)>1 || max(p)<0
       error('p needs to be between 0 and 1 (probability values)'); 
    end

    numCells = size(p,1);
    N = numCells*numCells-numCells;
    p(logical(eye(size(p)))) = nan; % ignore diagonal
    p_avg = nanmean(p(:));

    % Two Neuron Connectivity Patterns
    % p is probability of being connected
    %   unidirectionally connected: p*(1-p) + (1-p)*p
    %   bidirectionally connected: p^2
    %   unconnected: (1-p)^2
    % unidirectional, bidirectional, unconnected
    doubletMotifProbability_avg = [2*p_avg*(1-p_avg) p_avg^2 (1-p_avg)^2];

    %%
    twoMotifs = nan(N,3);
    pSampleMean = nan(N,2); 
    c = 1; 
    for i = 1:numCells

        for j = 1:numCells

            if i==j
               continue; 
            end

            p1 = p(i,j);
            p2 = p(j,i);

            twoMotifs(c,:) = [p1*(1-p2)+p2*(1-p1) p1*p2 (1-p1)*(1-p2)]; 
            pSampleMean(c,:) = [p1 p2];
            c = c+1; 
        end
        fprintf('%d / %d\n',i,numCells);
    end

    %%
    doubletMotifProbability = mean(twoMotifs);
    doubletMotifLabel = {'uni','bi','unconn'}; 
    dev = doubletMotifProbability./doubletMotifProbability_avg; 
    p_avg_sample = mean(pSampleMean(:)); 
    doubletMotifProbability_avgSample = ...
            [2*p_avg_sample*(1-p_avg_sample) p_avg_sample^2 (1-p_avg_sample)^2];
        
    result.doubletMotifProbability = doubletMotifProbability;
    result.doubletMotifLabel = doubletMotifLabel; 
    result.deviation = dev; 
    result.doubletMotifProbability_avg = doubletMotifProbability_avg;
    result.p_avg_sample = p_avg_sample; 
    result.doubletMotifProbability_avgSample = doubletMotifProbability_avgSample; 
end