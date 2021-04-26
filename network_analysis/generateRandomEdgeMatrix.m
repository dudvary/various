function edgeMatrix = generateRandomEdgeMatrix( ...
                                        numNodes,numDirEdges,numUndirEdges)
% Generates binary edge matrix having numEdge number of Edges (without
% counting edges on diagonal)
% edgeMatrix is of size numNodes and constraint by 
%   (1) numUndirEdges: number of undirected edges 
%   (2) numDirEdges: number of direcnumNoted edges
%   each edgeMatrix has many combinations of possible number of directed
%   edges that match number of undirected edges
% Input:
% - numNodes: number of nodes (size of square matrix)
% - numDirEdge: edgeMatrix with numDirEdges of size numNodes x numNodes 
% - numUndirEdges: edgeMatrix with numUndirEdges of size numNodes x numNodes   
% Output:
% - edgeMatrix: binary matrix [numNodes x numNodes]
%               with numUndirEdges and numDirEdge

    maxNumDirEdges = (numNodes*numNodes-numNodes);
    
    if numNodes<=1
       error('numNodes has to be larger than 1!'); 
    end
    
    if numDirEdges<0 || numUndirEdges<0
       error('Number of directed and undirected edges has to be a positive number!'); 
    end

    if numDirEdges>maxNumDirEdges
        error('Cannot create more edges (%d) than entries in matrix (%d)!', ...
            numDirEdges,numNodes*numNodes-numNodes);
    end
    
    if numUndirEdges>numDirEdges
        error(['Number of undirected edges (%d) cannot be larger than ' ...
            ' number of directed edges (%d)!'], ...
            numUndirEdges,numDirEdges);
    end
    
    % --------------------------------------------
    % SPECIAL CASE 1
    % fully connected
    if numDirEdges==maxNumDirEdges
        edgeMatrix = true(numNodes,numNodes);
        return;
    end
    % --------------------------------------------
    
    % --------------------------------------------
    % SPECIAL CASE 2
    % fully unconnected
    if numDirEdges==0 || numUndirEdges==0
        edgeMatrix = false(numNodes,numNodes);
        return;
    end    
    % --------------------------------------------
    
    % --------------------------------------------
    % ALL OTHER CASES
    % Extract indicies non-diagonal entries in matrix
    idDiagonal = 1:numNodes+1:(numNodes^2);
    idNonDiagonal = setdiff(1:(numNodes^2),idDiagonal); 
    
    if numel(idNonDiagonal)~=maxNumDirEdges
        error(['Something is wrong! Number of non-diagonal entries (%d)' ...
            ' does not match maximal number of directed edges (%d)!'], ...
            numel(idNonDiagonal),maxNumDirEdges);
    end
    
    % do not go through all possibilities at once, but draw each index
    % subsequently and then restrict the number of possibilities
    % accordingly
    % this is way faster and more efficient than simply randperm
    numFreeEdges = numDirEdges-numUndirEdges;
    edgeMatrix = false(numNodes,numNodes);
    idxValid = idNonDiagonal; 

    for i = 1:numDirEdges
        
        % Get index of one random node
        idx = randperm(length(idxValid),1);
        edgeMatrix(idxValid(idx)) = true; 
        
        % Find tranpose of drawn index
        [x,y] = ind2sub(size(edgeMatrix),idxValid(idx));
        idxForDel = sub2ind(size(edgeMatrix),y,x); 
        
        % check if transposed index is already 1
        % if so increase number of numFreeEdges by one
        % numFreeEdges is the number of edges that can be choosen fully
        % random, every edge after this can only be an edge / index whose
        % transpose is not yet taken
        if edgeMatrix(idxForDel) == 1
            numFreeEdges = numFreeEdges + 1; 
        end
        
        % Delete drawn index
        idxValid(idx) = [];
        
        if (i-1)>=numFreeEdges 
            % Delete its transpose
            idxValid(idxValid==idxForDel) = []; 
        end
    end

    % Double check whether it worked and the number of edges match
    t = tril(edgeMatrix | edgeMatrix',-1);
    currNumUndirEdges = sum(t(:)); 

    if currNumUndirEdges~=numUndirEdges
        error(['Failed generating proper edgeMatrix that matches the ' ... 
            'constraints [numNodes = %d, numDirEdges = %d, ' ...
            'numUndirEdges = %d]!'],numNodes,numDirEdges,numUndirEdges); 
    end
    
    if sum(edgeMatrix(:))~=numDirEdges
        error(['Failed generating proper edgeMatrix that matches the ' ... 
            'constraints [numNodes = %d, numDirEdges = %d, ' ...
            'numUndirEdges = %d]!'],numNodes,numDirEdges,numUndirEdges); 
    end
    % --------------------------------------------
end