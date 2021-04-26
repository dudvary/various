% [points] = readLandmarkFile(inputFilename)
% READ landmarkAscii file which can also be read by amira
% Input:
% - inputFilename: Location of .landmarkAscii
% Output:
% - points: returns points [ numPts x Dim ]
function [points] = readLandmarkFile(inputFilename)
        
    %% Open the text file.
    if (~(strcmp(inputFilename(end-13:end), '.landmarkAscii')))
        inputFilename = [inputFilename '.landmarkAscii'];
    end

    fileID = fopen(inputFilename,'r'); 

    if fileID==-1
       error(['Cannot read Amira File ' inputFilename]); 
    end

    %% Find Number of Points and there dimensionality (3)
    numPts = nan; 
    dim = nan; 
    tline = fgetl(fileID); 
    while ischar(tline) && (isnan(numPts) || isnan(dim))

        idx = strfind(tline,'define Markers'); 
        if ~isempty(idx)
            C = textscan(tline(idx:end),'%*s %*s %d','delimiter',' ');
            numPts = cell2mat(C);
        end

        idx = strfind(tline,'Markers { float['); 
        if ~isempty(idx)
            dim = str2double(tline(strfind(tline,'[')+1:strfind(tline,']')-1));
        end

        tline = fgetl(fileID); 
    end

    if isnan(numPts)
       error('No Definition of Markers were found! (Number of Points)'); 
    end

    if isnan(dim)
       error('No Definition of Dimensionality were found!'); 
    end

    %% Read in all points
    while ischar(tline) && strcmp(tline,'@1')==0
        tline = fgetl(fileID);    
    end

    if (strcmp(tline,'@1')==0)
       warning('No points could be found in landmarkAscii!');  
       points = [];
       return;
    end

    % numPts = numPts; 

    points = nan(numPts,dim); 
    i = 1; 
    while ischar(tline) && i < numPts+1

        tline = fgetl(fileID); 
        C = textscan(tline,'%f64');
        points(i,:) = cell2mat(C)';
        i = i + 1; 
    end

    tline = fgetl(fileID);

    while ischar(tline) && strcmp(tline,''); 
        tline = fgetl(fileID);
    end
    if ischar(tline)
       warning(['There is still input left: ' tline]) 
    end

    if sum(isnan(points(:)))>0 
       warning('There are still NaN points (i.e. some points have not been added!)'); 
    end
        

end
