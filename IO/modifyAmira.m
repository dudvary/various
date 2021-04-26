% [] = modifyAmira(inputfilename,outputfilename,scalingFactor)
% Reads in ImageDataVolume Amira file and writes the same Amira file, 
% but with scaled values
% Input:
% - inputfilename: /path/to/inputfilename.am
% - outputfilename: /path/to/outputfilename.am
% - scalingFactor: scaling Factor to multiply values with
function [] = modifyAmira(inputfilename,outputfilename,scalingFactor)

    if isempty(scalingFactor)
        error('Scaling factor is empty!')
    end

    % Open file.
    if (~(strcmp(inputfilename(end-2:end), '.am')))
        inputfilename = [inputfilename '.am'];
    end
    
    % Output file.
    if (~(strcmp(outputfilename(end-2:end), '.am')))
        outputfilename = [outputfilename '.am'];
    end
    
    fileID = fopen(inputfilename,'r');
    
    if fileID==-1
       error(['Cannot read Amira File ' inputfilename]); 
    end

    strTxt = cell(0); 
    tline = fgetl(fileID); 
    strTxt{end+1} = tline; 

    while ischar(tline) 
        tline = fgetl(fileID); 
        strTxt{end+1} = tline; 
        
        if strcmp(tline,'@1')
            break;
        end
    end
    
    while ischar(tline)
        tline = fgetl(fileID); 
        
        if strcmp(tline,'')
           break; 
        end
        
        C = textscan(tline,'%f64');
        tmp = scalingFactor.*cell2mat(C);
        strTxt{end+1} = num2str(tmp); 
    end
    
    % Close the text file.
    fclose(fileID);
    
    % Write Copy
    fid = fopen(outputfilename,'w');
    if fid ~= -1
        for i = 1:length(strTxt);            
            fprintf(fid,'%s\n',strTxt{i});
        end
        fclose(fid);
    else
       error(['Cannot write ' outputfilename '!']); 
    end

end