% [somaCenter] = findSomaCenterFromHoc(filename)
% Get 3D point coordinates of soma of SpatialGraph defined in .hoc file
% Input:
% - filename: /path/to/hocfile.hoc
% Output:
% - somaCenter: 3D point coordinates of soma
function [somaCenter] = findSomaCenterFromHoc(filename)

    % Open the text file.
    if (~(strcmp(filename(end-3:end), '.hoc')))
        filename = [filename '.hoc'];
    end
    
    fileID = fopen(filename,'r');
    
    tline = fgetl(fileID); 
    
    somaCenter = []; 
    
    while ischar(tline)
        
        idxSoma = strfind(tline,'create soma'); 
        
        % If soma definition found
        if (~isempty(idxSoma))
            
            while isempty(strfind(tline,'{pt3dadd('))
                tline = fgetl(fileID); 
            end
            
            % 3D Points have been found, store them
            while (~isempty(strfind(tline,'{pt3dadd(')))
                idxEnd = strfind(tline,')}'); 
                C = textscan(tline(10:idxEnd),'%f','delimiter',',');
                dim = cell2mat(C)';
                somaCenter = [somaCenter; dim(1:3)];
                tline = fgetl(fileID); 
            end
            
            break;             
        end
        
        tline = fgetl(fileID); 
    end

    if isempty(idxSoma)
       error('No Soma has been found!');  
    end
    
    fclose(fileID);
end