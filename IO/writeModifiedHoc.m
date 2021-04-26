% [] = writeModifiedHoc(treeDend,treeApic,filename,outputname)
% Reads in hoc file and modifies the diameters of this .hoc file 
% stores modified .hoc file under outputname
% Diameters are changed by matching the coordinates (X Y Z) and then simply
% changing the corresponding diameter. Mean over diameter is taking in case
% there are points with the same coordinates (starting and ending points of
% branching edges)
% Input:
% - treeDend: tree Structure of BasalDendrite
% - treeApic: tree Structure of ApicalDendrite
% - filename: path/to/file.hoc
% - outputname: path/to/outputhoc
function [] = writeModifiedHoc(treeDend,treeApic,filename,outputname)

    % Append .hoc if missing
    if (~(strcmp(filename(end-3:end), '.hoc')))
        filename = [filename '.hoc'];
    end
    
    % Append .hoc if missing
    if (~(strcmp(outputname(end-3:end), '.hoc')))
        outputname = [outputname '.hoc'];
    end

    % Read File
    fileID = fopen(filename,'r');
    if (fileID==-1)
       error(['Cannot open ' filename '!']); 
    end

    strTxt = cell(0); 
    newDia = []; 
    tline = fgetl(fileID); 
    c = 1; 
    SegFound = 0; 
    
    while ischar(tline)
        
        % Copy strings
        strTxt{end+1} = tline; 
        
        if SegFound == 0
            SegFound = ~isempty(strfind(tline,'dend')) || ~isempty(strfind(tline,'apical')) ; 
        end
        
        % Only check if dend or apical have been found (skip soma)
        if SegFound && ~isempty(strfind(tline,'{pt3dadd(')) && isempty(strfind(tline,'Spine'))

                    idxEnd = strfind(tline,')}'); 
                    C = textscan(tline(10:idxEnd),'%f64','delimiter',',');
                    pts = cell2mat(C)';

                    tmpDend = ismember(treeDend.X,pts(1)) + ...
                            ismember(treeDend.Y,pts(2)) + ...
                            ismember(treeDend.Z,pts(3));
                    
                    if ~isnan(treeApic)
                        tmpApic = ismember(treeApic.X,pts(1)) + ...
                            ismember(treeApic.Y,pts(2)) + ...
                            ismember(treeApic.Z,pts(3));
                    end
                    
                    % Takes the mean over all existing diameters for this
                    % point
                    if ~isnan(treeApic)
                        meantmp = mean([treeDend.D(tmpDend==3)' ... 
                            treeApic.D(tmpApic==3)']); 
                    else
                        meantmp = mean([treeDend.D(tmpDend==3)']); 
                    end
                        
                    newDia = [newDia; c meantmp]; 
        end
        
        tline = fgetl(fileID);  
        c = c+1; 
    end
    
    fclose(fileID);

    % Write Copy
    fid = fopen(outputname,'w');
    if fid ~= -1
        
        for i = 1:length(strTxt); 
            
            % If index has a modified diameter, extract new diameter
            if sum(ismember(newDia(:,1),i)==1)
                tmpDia = newDia(newDia(:,1)==i,2);
                
                % Write new diameter instead of old diameter
                idxStart = strfind(strTxt{i},',');
                strTxt{i} = [strTxt{i}(1:idxStart(end)) ' ' num2str(tmpDia)  ')}']; 
            end
            
            fprintf(fid,'%s\n',strTxt{i});
        end
        fclose(fid);
    else
       error(['Cannot write ' outputname '!']); 
    end
end
