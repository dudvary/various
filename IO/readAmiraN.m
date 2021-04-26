% [matrix, origin, VoxelSize] = readAmiraN(filename)
% Read simple ImageDataVolume Amira file with a Number of Scalar Components
% larger than 1.
% Returns a 4D-matrix and coorindates of origin
% Input:
% - filename: /path/to/filename
% Output:
% - matrix: [N1 x N2 x N3 x N4] 4D matrix containing values of ImageDataVolume
% - origin: [1 x 3] relativ position of the origin in N1 N2 N3
% - VoxelSize: [1 x 3] size of each voxel (typically 50 x 50 x 50)
function [matrix, origin, VoxelSize] = readAmiraN(filename)

    % Open the text file.
    if (~(strcmp(filename(end-2:end), '.am')))
        filename = [filename '.am'];
    end
    
    fileID = fopen(filename,'r');
    
    % Compute dimensionality and BoundingBox of input matrix
    dim = nan(1,3); 
    BoundingBox = nan(1,6); 
    numScalar = 1; 
    tline = fgetl(fileID); 
    flag = 0; 
    
    while (ischar(tline) && flag<3)
        
        idx = strfind(tline,'define Lattice'); 
        if ~isempty(idx)
            C = textscan(tline(idx:end),'%*s %*s %d %d %d','delimiter',' ');
            dim = cell2mat(C);
            flag = flag+1; 
        end
        
        idx = strfind(tline,'BoundingBox'); 
        if ~isempty(idx)
            C = textscan(tline(idx:end),'%*s %d %d %d %d %d %d','delimiter',' ');
            BoundingBox = cell2mat(C);
            flag = flag+1; 
        end
        
        idx = strfind(tline,'{ float');
        if ~isempty(idx)
            flag = flag+1; 
            if strfind(tline,'[') % Otherwise scalar is 1
                C = textscan(tline(idx:end),'%*s %d %*s','delimiter',{'[',']'});
                numScalar = double(cell2mat(C));
            end
        end
        tline = fgetl(fileID); 
    end
    
    if sum(isnan(dim)>0)
       error('Definition of Lattice could not be found!');  
    end
    
    if sum(isnan(BoundingBox)>0)
       error('Definition of BoundingBox could not be found!');  
    end
    
    BoundingBox = double(BoundingBox);   
    if numScalar>1
        matrix = zeros(dim(1),dim(2),dim(3),numScalar); 
    else
        matrix = zeros(dim); 
    end
    idx = 0; 
    
    % Fill matrix with the corresponding values
    while (ischar(tline) && idx <= numel(matrix))
        
        if strcmp(tline,'@1')
            idx = idx+1;    
            tline = fgetl(fileID); 
        end
        
        if (idx > 0 && idx <= numel(matrix))
            C = textscan(tline,'%f64');
            tmp = cell2mat(C); 
            if length(tmp)~=numScalar
                error(['Problem with number of scalar components! ' ...
                    num2str(length(tmp)) ' data points were found. Only ' ...
                    num2str(numScalar) ' have been defined!'])
            end
            
            if numScalar>1
                [x,y,z] = ind2sub(dim,idx); 
                matrix(x,y,z,:) = tmp; 
            else
                matrix(idx) = tmp; 
            end
            
            idx = idx + 1; 
        end
        tline = fgetl(fileID); 
    end
    
    if ischar(tline)
       warning('Something did not work! There are still elements in the list that should be inserted into the matrix.'); 
    end
    
    % Compute VoxelSize
    VoxelSize(1) = abs(BoundingBox(1)-BoundingBox(2))/(dim(1)-1); 
    VoxelSize(2) = abs(BoundingBox(3)-BoundingBox(4))/(dim(2)-1); 
    VoxelSize(3) = abs(BoundingBox(5)-BoundingBox(6))/(dim(3)-1); 
    
    if VoxelSize(1) ~= VoxelSize(2) || VoxelSize(2) ~= VoxelSize(3)
       warning(['VoxelSize is not the same for all dimensionalities! ' ...
           num2str(VoxelSize)]) 
    end
    VoxelSize = double(VoxelSize); 

    % Compute Origin
    origin = BoundingBox([1 3 5])./VoxelSize; 
    origin = abs(origin) + 1; 

    % Close the text file.
    fclose(fileID);
end