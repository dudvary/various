%% Returns motif and pMotif for each of the 16 motifs
% Input:
% - MotifSpectrum: Output from getMotifSpectrum();
% - idxMotifClass: [1 16]
% - p: from \data\Sampling_C2.mat NumTrails x 64
% - cell1: CellType of Cell1
% - cell2: CellType of Cell2
% - cell3: CellType of Cell3
% - preCellsList: List of PresynapticCellTypes (10 + VPM)
% Output
% - motif
%   - motif.name: CellTypes of cell1,cell2, and cell3
%   - motif.pType: avg. probability of CellType Combination having this
%               motif
%   - motif.p: sum of probabilities
%   - motif.n: sample size
%   - motif.str: connectivity of triplet motif
% - pMotif: normalized sum of probabilities (motif.p), sampled probability of
%               certain triplet motif occuring           
function [motif,pMotif] = motifSampling(MotifSpectrum,idxMotifClass, ...
                                        p,cell1,cell2,cell3,preCellsList)

    idxMotif = find(MotifSpectrum.Class==idxMotifClass);
    numSamples = size(p,1); 
    
    if isempty(idxMotif)
       error(['MotifClass ' num2str(idxMotifClass) ' was not found!']); 
    end
    
    pValues = nan(numSamples*length(idxMotif),1);
    
    for i = 1:length(idxMotif)
        pValues(1+(i-1)*numSamples:numSamples*i) = p(:,idxMotif(i)); 
    end
    
    switch idxMotifClass
        case 1 % 1 Motif
            [motif,pMotif] = motifSampling_OrderIrrelevant(pValues,cell1,cell2,cell3,preCellsList);
            motif.str = '1<->2<->3<->'; 
        case 2 % 2 Motifs
            [motif,pMotif] = motifSampling_2(pValues,cell1,cell2,cell3,preCellsList);
        case 3 % 6 Motifs
            [motif,pMotif] = motifSampling_3(pValues,cell1,cell2,cell3,preCellsList);
        case 4 % 6 Motifs
            [motif,pMotif] = motifSampling_4(pValues,cell1,cell2,cell3,preCellsList);
        case 5 % 3 Motifs
            [motif,pMotif] = motifSampling_5(pValues,cell1,cell2,cell3,preCellsList); 
        case 6 % 3 Motifs
            [motif,pMotif] = motifSampling_6(pValues,cell1,cell2,cell3,preCellsList); 
        case 7 % 6 Motifs
            [motif,pMotif] = motifSampling_7(pValues,cell1,cell2,cell3,preCellsList); 
        case 8 % 3 Motifs
            [motif,pMotif] = motifSampling_8(pValues,cell1,cell2,cell3,preCellsList); 
        case 9 % 6 Motifs
            [motif,pMotif] = motifSampling_9(pValues,cell1,cell2,cell3,preCellsList); 
        case 10 % 6 Motifs
            [motif,pMotif] = motifSampling_10(pValues,cell1,cell2,cell3,preCellsList); 
        case 11 % 6 Motifs
            [motif,pMotif] = motifSampling_11(pValues,cell1,cell2,cell3,preCellsList); 
        case 12 % 3 Motifs
            [motif,pMotif] = motifSampling_12(pValues,cell1,cell2,cell3,preCellsList);  
        case 13 % 3 Motifs
            [motif,pMotif] = motifSampling_13(pValues,cell1,cell2,cell3,preCellsList);           
        case 14 % 3 Motifs
            [motif,pMotif] = motifSampling_14(pValues,cell1,cell2,cell3,preCellsList);
        case 15 % 6 Motifs
            [motif,pMotif] = motifSampling_15(pValues,cell1,cell2,cell3,preCellsList);
        case 16 % 1 Motifs
            [motif,pMotif] = motifSampling_OrderIrrelevant(pValues,cell1,cell2,cell3,preCellsList);
            motif.str = '1 2 3';
        otherwise
            error('MotifClass between 1 and 16');
    end
end

function [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList)
% O = n*n*n

    n = length(preCellsList); % 11
    numPossibilites = n*n*n;
    
    motif.name = cell(numPossibilites,3);
    motif.n = nan(numPossibilites,1);
    motif.p = nan(numPossibilites,1); 
    c = 1; 

    for t1 = 1:length(preCellsList)

        idx1 = strcmp(cellArray(:,1),preCellsList{t1});
        
        if sum(idx1)==0
            fprintf(['WARNING! No Combination found for ' preCellsList{t1} '\n']); 
            continue;
        end

        for t2 = 1:length(preCellsList)
            
            idx2 = strcmp(cellArray(:,2),preCellsList{t2});
            
            if sum(idx1 & idx2)==0
                fprintf(['WARNING! No Combination found for ' preCellsList{t1} ...
                        ' ' preCellsList{t2} '\n']); 
                continue;
            end

            for t3 = 1:length(preCellsList)

                idx3 = strcmp(cellArray(:,3),preCellsList{t3});

                idxFinal = idx1 & idx2 & idx3;
                
                if (sum(idxFinal)==0)
                    %fprintf(['WARNING! No Combination found for ' preCellsList{t1} ...
                    %    ' ' preCellsList{t2} ' ' preCellsList{t3} '\n']); 
                    continue
                end
                
                nametmp = cellArray(idxFinal,:);
                motif.name(c,:) = nametmp(1,:);
                motif.n(c) = sum(idxFinal);
                motif.p(c) = sum(pValues(idxFinal));
                motif.pType(c) = mean(pValues(idxFinal));

                c = c+1; 
            end
        end
    end

    if (c-1 ~= numPossibilites)
       warning(['Only ' num2str(c-1) ' out of ' num2str(numPossibilites) ...
           ' possibilities checked!']); 
    end
    
    pMotif = motif.p./sum(motif.p);
end

function [numPossibilites] = loopComputation(n,k)

    if nargin<2
        k = 3;
    end

    if n<170 % values above are inf
        numPossibilites = factorial(k+n-1)/(factorial(k)*factorial(n-1));
    else
        numPossibilites = 0; 
        for t1 = 1:n 
            for t2 = 1:n 
                if (t1>t2) 
                    continue;
                end

                for t3 = 1:n 
                    if (t2>t3)
                        continue;
                    end
                    numPossibilites = numPossibilites+1;
                end
            end
        end
    end
end

function [motif,pMotif] = allCombinationsIrrelevantOrdering(cellArray,pValues,preCellsList)
%     
    n = length(preCellsList);
    numPossibilites = loopComputation(n,3);
    
    motif.name = cell(numPossibilites,3);
    motif.n = nan(numPossibilites,1);
    motif.p = nan(numPossibilites,1); 
    c = 1; 

    for t1 = 1:length(preCellsList)
        
        idx1_1 = strcmp(cellArray(:,1),preCellsList{t1});
        idx2_1 = strcmp(cellArray(:,2),preCellsList{t1});
        idx3_1 = strcmp(cellArray(:,3),preCellsList{t1});

        for t2 = 1:length(preCellsList)
            
            if (t1>t2)
                continue;
            end

            idx1_2 = strcmp(cellArray(:,1),preCellsList{t2});            
            idx2_2 = strcmp(cellArray(:,2),preCellsList{t2}); 
            idx3_2 = strcmp(cellArray(:,3),preCellsList{t2}); 
            
            idx1_1_2_2 = idx1_1 & idx2_2; 
            idx1_1_3_2 = idx1_1 & idx3_2;
            idx1_2_2_1 = idx1_2 & idx2_1;
            idx1_2_3_1 = idx1_2 & idx3_1; 
            idx2_1_3_2 = idx2_1 & idx3_2;
            idx2_2_3_1 = idx2_2 & idx3_1;

            for t3 = 1:length(preCellsList)

                if (t2>t3)
                    continue;
                end
                
                idx1 = idx1_1_2_2 & strcmp(cellArray(:,3),preCellsList{t3});
 
                idx2 = idx1_1_3_2 & strcmp(cellArray(:,2),preCellsList{t3});
                
                idx3 = idx1_2_2_1 & strcmp(cellArray(:,3),preCellsList{t3});
                    
                idx4 = idx1_2_3_1 & strcmp(cellArray(:,2),preCellsList{t3});
                    
                idx5 = idx2_1_3_2 & strcmp(cellArray(:,1),preCellsList{t3});
                    
                idx6 = idx2_2_3_1 & strcmp(cellArray(:,1),preCellsList{t3});
                    
                idxFinal = idx1 | idx2 | idx3 | idx4 | idx5 | idx6;
                
                if (sum(idxFinal)==0)
                    %fprintf(['WARNING! No Combination found for ' preCellsList{t1} ...
                    %    ' ' preCellsList{t2} ' ' preCellsList{t3} '\n']); 
                    continue
                end
                
                nametmp = [cellArray(idxFinal,1) cellArray(idxFinal,2) cellArray(idxFinal,3)];
                motif.name(c,:) = nametmp(1,:);
                motif.n(c) = sum(idxFinal);
                motif.p(c) = sum(pValues(idxFinal));
                motif.pType(c) = mean(pValues(idxFinal));

                c = c+1; 
            end
        end
        
        fprintf('%d/%d computed\n',t1,n);
        
    end

    if (c-1 ~= numPossibilites)
       warning(['Only ' num2str(c-1) ' out of ' num2str(numPossibilites) ...
           ' possibilities checked!']); 
    end
    
    pMotif = motif.p./sum(motif.p);
end

function [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList)
    % for cellArray(:,2) and cellArray(:,3) order is irrelevant
    % cellArray(:,1) is fixed
    
    % O = n * ((n*n-n)/2+n);
    n = length(preCellsList); % 11
    numPossibilites = n*((n*n-n)/2+n); % 726 (n^3+n^2)/2

    motif.name = cell(numPossibilites,3);
    motif.n = nan(numPossibilites,1);
    motif.p = nan(numPossibilites,1); 
    c = 1; 

    for fixedType = 1:length(preCellsList)

        idxFixed = strcmp(cellArray(:,1),preCellsList{fixedType});
        
        if sum(idxFixed)==0
            fprintf(['WARNING! No Combination found for ' preCellsList{fixedType} '\n']); 
            continue
        end

        for interchangeableType1 = 1:length(preCellsList)
            
            idx1_3 = strcmp(cellArray(:,3),preCellsList{interchangeableType1});
            idx1_2 = strcmp(cellArray(:,2),preCellsList{interchangeableType1});

            for interchangeableType2 = 1:length(preCellsList)

                if (interchangeableType1>interchangeableType2)
                    continue;
                end
                
                idx1 = idx1_3 & strcmp(cellArray(:,2),preCellsList{interchangeableType2});
                idx2 = idx1_2 & strcmp(cellArray(:,3),preCellsList{interchangeableType2});

                idxFinal = idxFixed & (idx1 | idx2);
                
                if (sum(idxFinal)==0)
                    %fprintf(['WARNING! No Combination found for ' preCellsList{fixedType} ...
                    %    ' ' preCellsList{interchangeableType1} ' ' preCellsList{interchangeableType1} '\n']); 
                    continue
                end
                
                nametmp = cellArray(idxFinal,:);
                motif.name(c,:) = nametmp(1,:);
                motif.n(c) = sum(idxFinal);
                motif.p(c) = sum(pValues(idxFinal));
                motif.pType(c) = mean(pValues(idxFinal));

                c = c+1; 
            end
        end
    end

    pMotif = motif.p./sum(motif.p);

end

% Motif 1 and motif 16
function [motif,pMotif] = motifSampling_OrderIrrelevant(pValues,cell1,cell2,cell3,preCellsList)
    % motif1:  1<->2<->3<->
    % motif16: 1 2 3
    cellArray = [cell1 cell2 cell3];       
    [motif,pMotif] = allCombinationsIrrelevantOrdering(cellArray,pValues,preCellsList);
end

function [motif,pMotif] = motifSampling_2(pValues,cell1,cell2,cell3,preCellsList)

    % 1->2->3->
    cellArray = [cell1 cell2 cell3];              % 1->2 2->3 3->1   1->2->3->
    cellArray = [cellArray; cell1 cell3 cell2];   % 1->3 2->1 3->2   1->3->2->  
    [motif,pMotif] = allCombinationsIrrelevantOrdering(cellArray,pValues,preCellsList);
    motif.str = '1->2->3->';    
end

function [motif,pMotif] = motifSampling_3(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1<->2<->3->  
    cellArray = [cell2 cell3 cell1];              % 1->2 1->3 2->3 3->1 3->2    2<->3<->1->
    cellArray = [cellArray; cell1 cell3 cell2];   % 1->3 2->1 2->3 3->1 3->2    1<->3<->2->
    cellArray = [cellArray; cell3 cell1 cell2];   % 1->2 1->3 2->1 2->3 3->1    3<->1<->2->
    cellArray = [cellArray; cell2 cell1 cell3];   % 1->2 1->3 2->1 3->1 3->2    2<->1<->3->
    cellArray = [cellArray; cell1 cell2 cell3];   % 1->2 2->1 2->3 3->1 3->2    1<->2<->3->
    cellArray = [cellArray; cell3 cell2 cell1];   % 1->2 1->3 2->1 2->3 3->2    3<->2<->1->
    
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1<->2<->3->';
end

function [motif,pMotif] = motifSampling_4(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1<->2->3->
    cellArray = [cell3 cell1 cell2];              % 1->2 1->3 2->3 3->1    3<->1->2->
    cellArray = [cellArray; cell1 cell3 cell2];   % 1->3 2->1 3->1 3->2    1<->3->2->
    cellArray = [cellArray; cell1 cell2 cell3];   % 1->2 2->1 2->3 3->1    1<->2->3->
    cellArray = [cellArray; cell2 cell1 cell3];   % 1->2 1->3 2->1 3->2    2<->1->3->
    cellArray = [cellArray; cell2 cell3 cell1];   % 1->2 2->3 3->1 3->2    2<->3->1->
    cellArray = [cellArray; cell3 cell2 cell1];   % 1->3 2->1 2->3 3->2    3<->2->1->
    
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1<->2->3->';
end

function [motif,pMotif] = motifSampling_5(pValues,cell1,cell2,cell3,preCellsList)

    % idx1 is fixed, idx2, and idx3 are interchangeable
    % 1<-2<->3->
    cellArray = [cell1 cell2 cell3];              % 2->1 2->3 3->1 3->2     1<-2<->3->           
    cellArray = [cellArray; cell2 cell3 cell1];   % 1->2 1->3 3->1 3->2     2<-3<->1->
    cellArray = [cellArray; cell3 cell1 cell2];   % 1->2 1->3 2->1 2->3     3<-1<->2->
    
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '1<-2<->3->';     
end

function [motif,pMotif] = motifSampling_6(pValues,cell1,cell2,cell3,preCellsList)

    % idx1 is fixed, idx2, and idx3 are interchangeable
    % <-1->2<->3
    cellArray = [cell1 cell2 cell3];              % 1->2 1->3 2->3 3->2     <-1->2<->3           
    cellArray = [cellArray; cell2 cell3 cell1];   % 1->3 2->1 2->3 3->1     <-2->3<->1
    cellArray = [cellArray; cell3 cell2 cell1];   % 1->2 2->1 3->1 3->2     <-3->2<->1
    
    % first is fixed, 2 and 3 are interchangeable
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '<-1->2<->3';   
end

function [motif,pMotif] = motifSampling_7(pValues,cell1,cell2,cell3,preCellsList)
    
    % <-1->2->3   
    cellArray = [cell1 cell2 cell3];              % 1->2 1->3 2->3    <-1->2->3
    cellArray = [cellArray; cell1 cell3 cell2];   % 1->2 1->3 3->2    <-1->3->2
    cellArray = [cellArray; cell2 cell1 cell3];   % 1->3 2->1 2->3    <-2->1->3
    cellArray = [cellArray; cell2 cell3 cell1];   % 2->1 2->3 3->1    <-2->3->1
    cellArray = [cellArray; cell3 cell1 cell2];   % 1->2 3->1 3->2    <-3->1->2
    cellArray = [cellArray; cell3 cell2 cell1];   % 2->1 3->1 3->2    <-3->2->1
    
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '<-1->2->3';   
end

function [motif,pMotif] = motifSampling_8(pValues,cell1,cell2,cell3,preCellsList)

    % 1<->2 3<->  
    cellArray = [cell1 cell2 cell3];              % 1->2 1->3 2->1 3->1   1<->2 3<->
    cellArray = [cellArray; cell3 cell2 cell1];   % 1->3 2->3 3->1 3->2   3<->2 1<-> 
    cellArray = [cellArray; cell2 cell3 cell1];   % 1->2 2->1 2->3 3->2   2<->3 1<->
        
    % first is fixed, 2 and 3 are interchangeable
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '1<->2 3<->';   
end

function [motif,pMotif] = motifSampling_9(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1->2->3   
    cellArray = [cell1 cell2 cell3];              % 1->2->3
    cellArray = [cellArray; cell3 cell2 cell1];   % 3->2->1
    cellArray = [cellArray; cell2 cell3 cell1];   % 2->3->1
    cellArray = [cellArray; cell1 cell3 cell2];   % 1->3->2
    cellArray = [cellArray; cell3 cell1 cell2];   % 3->1->2
    cellArray = [cellArray; cell2 cell1 cell3];   % 2->1->3
                    
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1->2->3';   
end

function [motif,pMotif] = motifSampling_10(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1<->2<-3   
    cellArray = [cell2 cell1 cell3];              % 2<->1->3
    cellArray = [cellArray; cell1 cell2 cell3];   % 1<->2->3
    cellArray = [cellArray; cell3 cell2 cell1];   % 3<->2->1
    cellArray = [cellArray; cell2 cell3 cell1];   % 2<->3->1
    cellArray = [cellArray; cell3 cell1 cell2];   % 3<->1->2
    cellArray = [cellArray; cell1 cell3 cell2];   % 1<->3->2
            
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1<->2->3'; 
end

function [motif,pMotif] = motifSampling_11(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1<->2<-3   
    cellArray = [cell2 cell1 cell3];              % 2<->1<-3
    cellArray = [cellArray; cell1 cell2 cell3];   % 1<->2<-3
    cellArray = [cellArray; cell3 cell2 cell1];   % 3<->2<-1
    cellArray = [cellArray; cell2 cell3 cell1];   % 2<->3<-1
    cellArray = [cellArray; cell3 cell1 cell2];   % 3<->1<-2
    cellArray = [cellArray; cell1 cell3 cell2];   % 1<->3<-2
    
    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1<->2<-3'; 
end

function [motif,pMotif] = motifSampling_12(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1<-2 3->
    cellArray = [cell1 cell2 cell3];            % 2->1<-3     1<-2 3->
    cellArray = [cellArray; cell2 cell1 cell3]; % 1->2<-3     2<-1 3->
    cellArray = [cellArray; cell3 cell2 cell1]; % 1->3<-2     3<-2 1->
    
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '1<-2 3->';   
end

function [motif,pMotif] = motifSampling_13(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1->2 3<-
    cellArray = [cell1 cell2 cell3];             % 1->2 1->3    1->2 3<-
    cellArray = [cellArray; cell2 cell1 cell3];  % 2->1 2->3    2->1 3<-
    cellArray = [cellArray; cell3 cell1 cell2];  % 3->1 3->2    3->1 2<-
    
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '1->2 3<-';   
end

function [motif,pMotif] = motifSampling_14(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1 2<->3
    cellArray = [cell3 cell1 cell2];                % 3 1<->2 
    cellArray = [cellArray; cell1 cell2 cell3];     % 1 2<->3 
    cellArray = [cellArray; cell2 cell1 cell3];     % 2 1<->3
    
    [motif,pMotif] = allCombinationsOneConstant(cellArray,pValues,preCellsList); 
    motif.str = '1 2<->3';   
end

function [motif,pMotif] = motifSampling_15(pValues,cell1,cell2,cell3,preCellsList)
    
    % 1->2 3
    cellArray = [cell1 cell2 cell3];                % 1->2 
    cellArray = [cellArray; cell2 cell1 cell3];     % 2->1 
    cellArray = [cellArray; cell2 cell3 cell1];     % 2->3
    cellArray = [cellArray; cell3 cell2 cell1];     % 3->2
    cellArray = [cellArray; cell1 cell3 cell2];     % 1->3
    cellArray = [cellArray; cell3 cell1 cell2];     % 3->1

    [motif,pMotif] = allCombinations(cellArray,pValues,preCellsList);
    motif.str = '1->2 3';
end