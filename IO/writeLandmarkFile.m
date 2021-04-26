% [] = writeLandmarkFile(points,outputFilename)
% Write a landmark file which can be read by amira
% Input:
% - points: cloud of points numPtsx3
% - outputFilename: /path/to/outputfilename.landmarkAscii
function [] = writeLandmarkFile(points,outputFilename)
        
        fname = [outputFilename '.landmarkAscii'];
        fid = fopen(fname,'w'); 
        
        numPts = size(points,1); 
        
        if size(points,2) ~= 3
           error('Points are not 3D-coordinates!') 
        end
        
        str = sprintf(['# AmiraMesh 3D ASCII 2.0\n\n' ...
            'define Markers %d\n\n' ...
            'Parameters {\n' ...
            '\tNumSets 1,\n' ...
            '\tContentType "LandmarkSet"\n}\n\n' ...
            'Markers { float[3] Coordinates } @1\n\n' ...
            '# Data section follows\n' ...
            '@1'],numPts); 
        
        if fid ~= -1
            fprintf(fid,'%s\n',str);   % no \r
            fclose(fid);
        end
        
        dlmwrite(fname,points,'-append','delimiter',' ','newline','pc'); 
end
