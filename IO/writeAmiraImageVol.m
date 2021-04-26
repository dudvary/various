% [] = writeAmiraImageVol(matrix,outputFilename,origin,VoxelSize)
% Writes ImageDataVolume so that it can be read by Amira
% Input:
% - matrix: 3D matrix 
% - outputFilename: path/to/outputfilename.am
% - origin: [1 x 3] vector with 3 elements of origin in indices 
%       (default: [1 1 1])
% - VoxelSize: [1 x 3] Size of each voxel
%       (default: [50 50 50]) 
function [] = writeAmiraImageVol(matrix,outputFilename,origin,VoxelSize)

    if nargin<3
        origin = ones(1,3);
    end
    if nargin<4
        VoxelSize = [50 50 50];
    end

    origin = double(origin); 
    BoundingBox = nan(1,6);
    BoundingBox([1 3 5]) = (origin-1).*(-VoxelSize);  
    BoundingBox([2 4 6]) = BoundingBox([1 3 5])+(size(matrix)-1).*VoxelSize; 

    if (strcmp(outputFilename(end-2:end), '.am'))
        fname = outputFilename;  
    else
        fname = [outputFilename '.am'];
    end
    
    fid = fopen(fname,'w');

    % Create Header
    strHeader = sprintf(['# AmiraMesh 3D ASCII 2.0\n\n' ...
                'define Lattice %d %d %d\n\n' ...                    
                'Parameters {\n' ...
                '    Content "%dx%dx%d float, uniform coordinates",\n' ...
                '    BoundingBox %d %d %d %d %d %d\n' ...
                '    CoordType "uniform"\n' ...
                '}\n\n' ...
                'Lattice { float Data } @1\n\n' ...
                '# Data section follows\n' ...
                '@1'],size(matrix,1),size(matrix,2),size(matrix,3), ...
                size(matrix,1),size(matrix,2),size(matrix,3), ...
                BoundingBox(1),BoundingBox(2),BoundingBox(3), ...
                BoundingBox(4),BoundingBox(5),BoundingBox(6));

    % Write to file
    if fid ~= -1
        fprintf(fid,'%s\n',strHeader);
        % Values
        for idx = 1:numel(matrix)
            fprintf(fid,'%e\n',matrix(idx));           
        end        
        fclose(fid);
    else
        error(['Failed writing ' fname ' to the disk!']); 
    end
end