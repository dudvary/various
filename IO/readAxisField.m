function [AxisFieldVectors,BoundingBox,Dimensions] = readAxisField(filename)
% [AxisFieldVectors,BoundingBox,Dimensions] = readAxisField(filename)
% Read in AxisField and return AxisFieldVectors, BoundingBox and Dimensions/Size
% Input:
% - filename
% Output:
% - AxisFieldVectors [dim_x*dim_y*dim_z x 3]: contains axis field vectors
%       for each position in grid 
% - BoundingBox [1x6] [xmin xmax ymin ymax zmin zmax]: BoundingBox of grid
% - Dimensions [1x3] [dim_x dim_y dim_z]: Size of grid
% 
% Voxel size can be computed from BoundingBox and Dimensions
%   voxelsize = (BoundingBox(2)-BoundingBox(1))/(dim(1)) + 1;
%   voxelsize = (BoundingBox(4)-BoundingBox(3))/(dim(2)) + 1;
%   voxelsize = (BoundingBox(5)-BoundingBox(5))/(dim(3)) + 1;

    fileID = fopen(filename,'r');

    if fileID==-1
       error(['Cannot read Axis Vector Field: ' filename]); 
    end

    foundLattice = 0;
    foundBB = 0; 
    foundData = 0; 
    counter = 1; 
    tline = fgetl(fileID); 

    while ischar(tline)

        if ~foundLattice
            idx = strfind(tline,'define Lattice'); 
            if ~isempty(idx)
                C = textscan(tline(idx:end),'%*s %*s %d %d %d','delimiter',' ');
                Dimensions = cell2mat(C);
                AxisFieldVectors = nan(prod(Dimensions),3);
                foundLattice = 1; 
            end
        end

        if ~foundBB
            idx = strfind(tline,'BoundingBox'); 
            if ~isempty(idx)
                C = textscan(tline(idx:end),'%*s %f %f %f %f %f %f','delimiter',' ');
                BoundingBox = cell2mat(C);
                foundBB  = 1; 
            end
        end

        if ~foundData && foundBB && foundLattice 
            idx = strcmp(tline,'@1'); 
            if idx==1
                foundData=1;
                tline = fgetl(fileID);
            end
        end

        if foundData
            C = textscan(tline,'%f %f %f','delimiter',' ');
            AxisFieldVectors(counter,:) = cell2mat(C); 
            
%             if sum(isnan(AxisFieldVectors(counter,:)))>0
%                 disp(tline);
%             end
            
            counter = counter+1; 
        end

        tline = fgetl(fileID); 
    end
    fclose(fileID);

    if foundLattice==0 || foundBB==0 || foundData==0
       error('Something went wrong! Lattice, BoundingBox or DataSection not found!'); 
    end
    
end