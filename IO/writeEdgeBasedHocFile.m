function [] = writeEdgeBasedHocFile(filename,pts)
% [] = writeEdgeBasedHocFile(filename,pts)
% Writes each cell in cell array pts as one edge to .hoc file
% Hoc file does not have any connectivity!
% Input:
% - filename: path/to/outputfile.hoc
% - pts: [1 x N] cell array
%       - each cell: [N x 4] 3D points array plus diameter value
% Output:
% - stores the EdgeBasedHocFile in filename

    % Append .hoc if missing
    if (~(strcmp(filename(end-3:end), '.hoc')))
        filename = [filename '.hoc'];
    end
    
    fileID = fopen(filename,'w');

    if (fileID==-1)
       error(['Cannot write to ' filename '!']); 
    end
    
    % Write Header (Comments)
    fprintf(fileID,'/*--------------------------------------------*/\n');
    fprintf(fileID,'/* Edge-based SpatialGraph (no connectivity!) */\n');
    fprintf(fileID,'/* created with writeEdgeBasedHocFile.m       */\n');
    fprintf(fileID,'/*--------------------------------------------*/\n');
    fprintf(fileID,'\n');
    
    % Write Edges
    for i = 1:length(pts)
        fprintf(fileID,'{create axon_%d}\n',i-1);
        fprintf(fileID,'{access axon_%d}\n',i-1);
        fprintf(fileID,'{nseg = 1}\n');
        fprintf(fileID,'{strdef color color = "Blue"}\n');
        fprintf(fileID,'{pt3dclear()}\n'); 
        fprintf(fileID,'/* Tree */\n'); 
        
        for j = 1:size(pts{i},1)
            tmpPt = pts{i}; 
            fprintf(fileID,'{pt3dadd(%f,%f,%f,%f)}\n',tmpPt(j,:));
        end
        fprintf(fileID,'\n');
    end
    
    fclose(fileID); 
end