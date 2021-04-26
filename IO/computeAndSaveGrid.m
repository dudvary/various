function [] = computeAndSaveGrid(BBmin,BBmax,filepath,voxelSize)
% computeAndSaveGrid(BBmin,BBmax,filepath,voxelSize)
% Write Grid
% Input:
% - BBmin: [1 x 3] minimal 3D Coordinates of Box
% - BBmax: [1 x 3] maximal 3D Coordinates of Box
% - filepath: path/to/outputfile.am
% - voxelSize: (optional, default: 50)

    if nargin<4
        voxelSize = 50; 
    end
    
    xValues = [BBmin(1):voxelSize:BBmax(1)];
    yValues = [BBmin(2):voxelSize:BBmax(2)];
    zValues = [BBmin(3):voxelSize:BBmax(3)];
    
    opt = 'xyz';

    for k = 1:3

        outputFilename = [filepath opt(k) '_' num2str(BBmin(k)) '_' ...
                            num2str(BBmax(k)) '.am'];
        pt = cell(0);

        switch k
            case 1 
                for y = 1:length(yValues)
                    for z = 1:length(zValues)
                        pt{end+1} = [xValues(1) yValues(y) zValues(z); ...
                                xValues(end) yValues(y) zValues(z)];
                    end
                end
            case 2
                for x = 1:length(xValues)
                    for z = 1:length(zValues)
                        pt{end+1} = [xValues(x) yValues(1) zValues(z); ...
                                xValues(x) yValues(end) zValues(z)];
                    end
                end
            case 3
                for y = 1:length(yValues)
                    for x = 1:length(xValues)
                        pt{end+1} = [xValues(x) yValues(y) zValues(1); ...
                                xValues(x) yValues(y) zValues(end)];
                    end
                end
        end

        writeSpatialGraphContour(pt,outputFilename);
    end
end