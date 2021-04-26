%% Extract Edge Statistics for given matrix
% Input:
% - pSelect: [numTrials x numCells x numCells] p values for all cell pairs
% - p_avg: average connectivity between all cells
% Output:
% - edgeStats
%     edgeStats.numCells
%     edgeStats.numOfTrials
%     edgeStats.pEdgeUniformMotifUndirEdges
%     edgeStats.pEdgeUniformMotifDirEdges
%     edgeStats.pEdgeMotifUndirEdges
%     edgeStats.pEdgeMotifDirEdges
%     edgeStats.freqOfDirEdges
%     edgeStats.freqOfUndirEdges  
%     edgeStats.listOfDirEdges
%     edgeStats.listOfUndirEdges
%     edgeStats.p_avg
function [edgeStats] = extractEdgeStatistics(pSelect,p_avg)

    numOfTrials = size(pSelect,1); 
    numCells = size(pSelect,2); 
    dispThreshold = 1e4; 
    
    % Square Matrix
    sz = [numCells numCells]; %size(p); 
    numDirEdges = sz(1)*sz(2)-min(sz); % number of directed edges
    numUndirEdges = numDirEdges/2; 
    
    fprintf('Matrix size = [%d %d]\n',sz(1),sz(2));
    fprintf('Number of Trials = %d\n',numOfTrials);
    fprintf('Number of directed Edges = %d\n',numDirEdges);
    fprintf('Number of undirected Edges = %d\n',numUndirEdges);
    % Matrix Mapper
    id = 1:prod(sz); % all ids
    idDiagonal = 1:sz(1)+1:prod(sz); % id of diagonal values
    idValid = setdiff(id,idDiagonal);
    if (numel(idValid)~=numDirEdges) % double check
        error('something is off');
    end

    listOfUndirEdges = 0:numUndirEdges;
    freqOfUndirEdges = zeros(size(listOfUndirEdges)); 
    pEdgeMotifUndirEdges = zeros(numOfTrials,size(listOfUndirEdges,2)); 
    pEdgeUniformMotifUndirEdges = zeros(size(listOfUndirEdges)); 

    listOfDirEdges = 0:numDirEdges;
    freqOfDirEdges = zeros(size(listOfDirEdges)); 
    pEdgeMotifDirEdges = zeros(numOfTrials,size(listOfDirEdges,2)); 
    pEdgeUniformMotifDirEdges = zeros(size(listOfDirEdges)); 
    
    prevProgress = 0; 

    % Go through all number of edges and combinations
    for currNumDirEdges = 0:numDirEdges

        if currNumDirEdges>numDirEdges
            error('%d has to be equal or smaller than %d', ...
                                            currNumDirEdges,numDirEdges);
        end

        ii = nchoosek(1:numDirEdges,currNumDirEdges); % indices with edges
        numCombinations = size(ii,1); % k

        fprintf('currNumEdges = %d : #Combinations = %d\n', ...
                                        currNumDirEdges,numCombinations);

        if numCombinations>dispThreshold
            fprintf('>> 0%%\n');
        end
                                    
        % index position in matrix
        idMatrix = idValid(ii);

        % now find number of undirected edges (use transpose feature)
        % think of a more efficient implementation!
        % Matrix Version
        for i = 1:numCombinations

            % Find undirected edges (1's in m or m')
            m = false(sz); 

            if ~isempty(idMatrix)
                if size(idMatrix,1)==1
                    if numCombinations==1 
                        % here all indices for one combination 
                        % (fully connected)
                        m(idMatrix) = 1; 
                    else        
                        % here each indices for one combination 
                        % (fully unconnected)
                        m(idMatrix(i)) = 1; 
                    end
                else
                    m(idMatrix(i,:)) = 1; 
                end
            end

            % Display progress
            if numCombinations>dispThreshold
                
                progress = floor(i/numCombinations*100);
                
                if (progress~=prevProgress || i==numCombinations)
                    if (floor(i/numCombinations*100))>9
                        str = '\b\b\b';
                    else
                        str = '\b\b';
                    end
                    fprintf([str '%d%%'],progress);  
                end
                
                prevProgress = progress;
            end
            
            % probability of combination
            pMotif = nan(1,numOfTrials); 
            for j = 1:numOfTrials
                p = squeeze(pSelect(j,:,:)); 
                pMotif(j) = prod(p(m & ~isnan(p))) * prod(1-p(~m & ~isnan(p)));
            end

            % probabiltiy of combination when average p
            connected = sum(m(:) & ~isnan(p(:)));
            unconnected = sum(~m(:) & ~isnan(p(:)));
            pMotifUniform = p_avg^connected * (1-p_avg)^unconnected;

            % ***************
            % Undirected Edges!
            % ***************
            % Number of Undirected Edges
            t = tril(m | m',-1);
            currNumUndirEdges = sum(t(:)); 
            idxUndirEdge = listOfUndirEdges==currNumUndirEdges; 
            freqOfUndirEdges(idxUndirEdge) = freqOfUndirEdges(idxUndirEdge) + 1; 

            % Motif Probability
            % connected x non-connected
            pEdgeMotifUndirEdges(:,idxUndirEdge) = ... 
                                pEdgeMotifUndirEdges(:,idxUndirEdge) + pMotif';

            % Motif Probability when average p only
            % connected x non-connected
            pEdgeUniformMotifUndirEdges(idxUndirEdge) = ...
                    pEdgeUniformMotifUndirEdges(idxUndirEdge) + pMotifUniform;

            % ***************
            % Directed Edges!
            % ***************
            idxDirEdge = listOfDirEdges==currNumDirEdges; 
            freqOfDirEdges(idxDirEdge) = freqOfDirEdges(idxDirEdge) + 1; 

            % Motif Probability
            % connected x non-connected
            pEdgeMotifDirEdges(:,idxDirEdge) = pEdgeMotifDirEdges(:,idxDirEdge) ...
                                            + pMotif';

            % Motif Probability when average p only
            % connected x non-connected
            pEdgeUniformMotifDirEdges(idxDirEdge) = ...
                        pEdgeUniformMotifDirEdges(idxDirEdge) + pMotifUniform;
        end
        
        if numCombinations>dispThreshold
            fprintf('\n');
        end
    end

    if abs(sum(pEdgeUniformMotifUndirEdges)-1)>1e-4
        error('Probability does not sum up to 1! (p = %.4f)', ...
                                sum(pEdgeUniformMotifUndirEdges));
    end

    % Create one struct for returning
    edgeStats.numCells = numCells; 
    edgeStats.pEdgeUniformMotifUndirEdges = pEdgeUniformMotifUndirEdges; 
    edgeStats.pEdgeUniformMotifDirEdges = pEdgeUniformMotifDirEdges; 
    edgeStats.pEdgeMotifUndirEdges = pEdgeMotifUndirEdges; 
    edgeStats.pEdgeMotifDirEdges = pEdgeMotifDirEdges; 
    edgeStats.freqOfDirEdges = freqOfDirEdges;
    edgeStats.freqOfUndirEdges = freqOfUndirEdges;    
    edgeStats.numOfTrials = numOfTrials;    
    edgeStats.listOfDirEdges = listOfDirEdges;  
    edgeStats.listOfUndirEdges = listOfUndirEdges;  
    edgeStats.p_avg = p_avg; 
end