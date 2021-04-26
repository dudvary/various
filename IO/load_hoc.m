% [tree] = load_hoc(filename,label,prefix_neurite)
% loads a .hoc file representing a tree and converts it into tree structure
% Tree structure as of https://www.treestoolbox.org/
% Cuntz H, Forstner F, Borst A, Häusser M (2010) 
% One Rule to Grow Them All: A General Theory of Neuronal Branching and 
% Its Practical Application. PLOS Computational Biology 6(8): e1000877. 
% https://doi.org/10.1371/journal.pcbi.1000877
% 
% Input 
% - filename: location and name of .hoc file
% - label (optional): [3 6 7 4] for [dend axon soma apical] 
% - neurite (optional): sometimes there's a prefix like 'marcel_dend', adds
% underscore to prefix automatically, set neurite to 'auto' and it detects
% prefix automatically by finding the soma and using the prefix
% Output
% - tree: tree cell array of size 4 for axon, dend, soma and apical. returns
%         empty cell if no tree has been found.  
% Improvements:
% - Fix Region Labeling, to distinguish between apical and basal dendrite
% - Fix prefix, just check for key words in line
% - Fix prefix, set to 'auto' to skip detect soma,dend,axon automatically
% - Fix apical tree output
% - Fix added warning if adjacency matrix dA is not square
% - Fix no warning when apical dendrite is missing
function [tree] = load_hoc(filename,label,prefix_neurite)

    % If only specific part of hoc should be load, like axon, dendrite, or
    % soma
    if nargin<2
        ward_list = [3 1 2]; % always start with soma to get soma coordinates
    else
        label(label==6) = 1; % axon
        label(label==3) = 2; % dend
        
        if sum(label>2)>0 || sum(label<1)>0
           error('Wrong label input! Use [3 6] for [dend axon soma]'); 
        end
        
        ward_list = [3 label];
    end
    
    if nargin<3
       prefix_neurite = '';  
    end
    
    % If Dendrite is defined, apical dendrit will be also extracted
    if ismember(2,ward_list) && ~ismember(4,ward_list)
       ward_list = [ward_list 4];  
    end

    % Append .hoc if missing
    if (~(strcmp(filename(end-3:end), '.hoc')))
        filename = [filename '.hoc'];
    end

    % Check whether prefix should be auto detected
    if strcmp(prefix_neurite,'auto')
        
        % Open File
        fileID = fopen(filename,'r');
        
        if (fileID==-1)
           error(['Cannot open ' filename '!']); 
        end
        
        tline = fgetl(fileID); 
        somaFound = 1; 
        
        while ischar(tline) && somaFound
            
            tmpStr = regexp(tline,'create \w*soma','match'); 
            
            if ~isempty(tmpStr) 
                if length(tmpStr{1})==11 % 'create soma' -> no prefix
                    prefix_neurite = ''; 
                else
                    prefix_neurite = tmpStr{1}(8:end-4); 
                end
                somaFound = 0; 
            end
            
            tline = fgetl(fileID); 
        end

        fclose(fileID);
        
        if somaFound==1
           error('No Soma was found. Could not identify prefix!') 
        end
    else % Adds Underscore to prefix_neurite
       prefix_neurite = [prefix_neurite '_'];  
    end
    
    % Create tree cell array (1: axon, 2: dendrite, 3: soma, 4: apical)
    tree = cell(1,4); 
    
    % Center of Soma (this is the root)
    somaCenter = nan(1,4); 
    
    % Tree Name
    treeName = strsplit(filename,{'/','\\'});
    
    % Store Soma, Axon, and Dendrite in seperate trees
    for ward = ward_list
        
        % Initialize tree
        tree{ward}.X = [];
        tree{ward}.Y = [];
        tree{ward}.Z = [];
        tree{ward}.dA = sparse (0); 
        tree{ward}.D = [];
        
        neurite_tmp = []; % Alternative neurite naming
        
        % Label arbors correctly
        switch ward
            case 1
                neuriteID = 6;
                neurite = 'axon';  
            case 2
                neuriteID = 3;
                neurite = 'dend'; 
                neurite_tmp = 'BasalDendrite';
            case 3
                neuriteID = 7;
                neurite = 'soma'; 
            case 4
                neuriteID = 4;
                neurite = 'apical'; 
        end
        
        if (nargin==3)
           neurite = [prefix_neurite neurite];  
        end
        
        tree{ward}.rnames = {num2str(neuriteID) neurite};
        tree{ward}.name = neurite; 
        
        % Open File
        fileID = fopen(filename,'r');
        
        if (fileID==-1)
           error(['Cannot open ' filename '!']); 
        end
        
        tline = fgetl(fileID); 
        
        % Connectivity:             tree{ward}.dA(2:end,2:end) = tmp; 

        % {1,:} : current name of edge
        % {2,:} : father point index
        % {3,:} : current point index
        % {4,:} : name of father edge
        Connectivity = cell(0);  
        ptIDX = 1; 
        
        while ischar(tline)

            % Check whether paragraph contains Axon, Dendrite, or Soma Segments
            SegFound = strfind(tline,['create ' neurite]); 
            
            if (isempty(SegFound)) && ischar(neurite_tmp)
                SegFound = strfind(tline,['create ' neurite_tmp]); 
            end

            % If Segment is found 
            if (~isempty(SegFound))

                if(neuriteID~=7) % Store Connectivity, name and pointID of first point
                    Connectivity{1,end+1} = tline(9:end-1); % name of edge
                    Connectivity{2,end} = ptIDX;            % point index of father edge
                end

                % Go until point definition
                while isempty(strfind(tline,'{pt3dadd('))
                    
                    ConFound = strfind(tline,['connect ' neurite]);
                    
                    if (isempty(ConFound)) && ischar(neurite_tmp)
                        ConFound = strfind(tline,['connect ' neurite_tmp]);
                    end
                    
                    if (~isempty(ConFound))
                        tmp = tline(strfind(tline,'), ')+3:strfind(tline,')}')-3);
                        tmpidx = strfind(tmp,'('); % Added because sometimes there's a double value in brackets
                        if ~isempty(tmpidx)
                            tmp = tmp(1:tmpidx(1)-1); 
                        end
                        Connectivity{4,end} = tmp; % name of father edge
                    end

                    tline = fgetl(fileID); 
                end

                % 3D Points have been found, store them. Last value is
                % diameter
                % All points in one segment are continous points
                segIDX = 1; 
                while (~isempty(strfind(tline,'{pt3dadd(')))

                    % Skip if spines are defined
                    while (~isempty(strfind(tline,'Spine')))
                        tline = fgetl(fileID); 
                    end
                    
                    % If Points are not defined in this line, segment is
                    % done. Happens if Spine is found as last point in
                    % segment.
                    if (isempty(strfind(tline,'{pt3dadd(')))
                        break;
                    end

                    idxEnd = strfind(tline,')}'); 
                    C = textscan(tline(10:idxEnd),'%f64','delimiter',',');
                    pts = cell2mat(C)';

                    % Add point coordinates
                    tree{ward}.X = [tree{ward}.X; pts(1)]; 
                    tree{ward}.Y = [tree{ward}.Y; pts(2)]; 
                    tree{ward}.Z = [tree{ward}.Z; pts(3)]; 

                    % Add diameter
                    tree{ward}.D = [tree{ward}.D; pts(4)]; 

                    if segIDX>1    % Skip first point, remaining points of 
                                   % one segment are all continoues points!
                        tree{ward}.dA(ptIDX,ptIDX-1) = 1; 
                    end

                    % Overall point count
                    ptIDX = ptIDX + 1;
                    
                    % Point count within segment
                    segIDX = segIDX + 1; 

                    tline = fgetl(fileID); 
                end

                % Store pointID of last point
                if(neuriteID~=7) 
                    Connectivity{3,end} = ptIDX-1; % Current point index
                else % Soma is only one connected segment, do not go through all other points
                    break; 
                end
                
            end

            tline = fgetl(fileID); 
        end

        fclose(fileID);

        % If soma, store center of soma plus mean diameter to use it as root
        if (neuriteID==7 && ~(isempty(tree{ward}.X))) 
            somaCenter = mean([tree{ward}.X tree{ward}.Y tree{ward}.Z tree{ward}.D]); 
        elseif sum(isnan(somaCenter))>0
            warning(['No soma has been found! Neurite will not be ' ...
                ' connected to soma/root!']);
        end
        
        % Check whether points have been found (sometimes dend is missing!)
        % If that's the case, report and skip to next tree
        if (isempty(tree{ward}.X) || isempty(tree{ward}.Y) || isempty(tree{ward}.Z))
            
            if neuriteID ~= 4
                warning(['No ' neurite ' coordinates have been found! ' ...
                    neurite ' is empty.'])
            end
            
            tree{ward} = []; 
            continue; 
        elseif (size(Connectivity,1) ~= 4 && neuriteID~=7)
            warning(['Something is wrong with the connectivity for ' neurite '!']); 
            tree{ward} = [];
            continue; 
        end
        
        % Boolean if soma was added (increases pts by 1) 
        addedPoint = sum(isnan(somaCenter))==0; 
        
        % Add Root point
        if (neuriteID~=7 && addedPoint)
            tree{ward}.X = [somaCenter(1); tree{ward}.X(:)];
            tree{ward}.Y = [somaCenter(2); tree{ward}.Y(:)];
            tree{ward}.Z = [somaCenter(3); tree{ward}.Z(:)];
            tree{ward}.D = [somaCenter(4); tree{ward}.D(:)];
            
            % Add last column of zeros (so that it is symmetric)
            % all points are made continous points in this step
            tmp = [tree{ward}.dA sparse(size(tree{ward}.dA,1),1)];
            tree{ward}.dA = sparse(size(tmp,1)+1,size(tmp,2)+1);
            tree{ward}.dA(2:end,2:end) = tmp; 
        else
            % Add last column of zeros (so that it is symmetric)
            % all points are made continous points in this step
            tree{ward}.dA = [tree{ward}.dA sparse(size(tree{ward}.dA,1),1)]; 
        end

        % Add Region Label
        tree{ward}.R = neuriteID.*ones(size(tree{ward}.D));
        
        % Get Adjacency Matrix correct
        % Add Branching Point connectivity
        if (neuriteID~=7) % Skip connectivity for for soma
                            
            for i = 1:size(Connectivity,2)

                currentName = Connectivity{1,i}; 
                currentPtIDX = Connectivity{3,i}+addedPoint; % Add because of root point

                for j = 1:size(Connectivity,2)
                    if strcmp(currentName,Connectivity{4,j}) 
                        tree{ward}.dA(Connectivity{2,j}+addedPoint,currentPtIDX) = 1; 
                    end
                end
                
                % If father is soma, establish connectivity to between
                % father point index and root idx 
                if (~isempty(strfind(Connectivity{4,i},'soma')) && ...
                        addedPoint)
                    tree{ward}.dA(Connectivity{2,i}+addedPoint,1) = 1; 
                end
                
            end
        end
        
        tree{ward}.ID = treeName{end}; 
        
        if size(tree{ward}.dA,1)~=size(tree{ward}.dA,2)
            warning('Adjacency Matrix is not square! Add one column of zeros ... Check this!'); 
            tree{ward}.dA = [tree{ward}.dA sparse(size(tree{ward}.dA,1),1)]; 
        end
        
    end  
end