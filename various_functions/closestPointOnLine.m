function [ptOnLine,tscale,distPtLine] = closestPointOnLine(pt,v1,v2)
% Compute the closest point on a line/vector  from another point
% Input:
% - pt: 2D or 3D point vector
% - v1: 2D or 3D point vector of line
% - v2: 2D or 3D point vector of line
% v1 and v2 define the line
% Output:
% - ptOnLine: pointOnLine
% - tscale: magnitude (relative positon on v1-v2) it's between 0 and 1
% - distPtLine: distance of point pt to line vector v1 and v2

    % If 2D, make it 3D by adding 0
    if (length(v1)==2 && length(v2)==2 && length(pt)==2)
       v1 = [v1 0];  
       v2 = [v2 0];  
       pt = [pt 0];      
    end

    if (length(v1)~=3 || length(v2)~=3 || length(pt)~=3)
        error('Vectors need to be same size. Only works for size 2 or 3!');
    end
    
    % all need to be row vectors!
    if iscolumn(v1)
        v1 = v1';
    end
    if iscolumn(v2)
        v2 = v2';
    end
    if iscolumn(pt)
        pt = pt';
    end   
    
    tscale = ((v1-v2)*(v2-pt)') / ((v1-v2) * (v1-v2)');
    tscale = -tscale;
    ptOnLine = v2 + tscale.*(v1-v2);
    
    if nargout>2
        distPtLine = sqrt(sum((ptOnLine - pt).^2));
        
        % For double checking
        distPtLine2 = distPointToLine(pt, v1, v2);
        
        if abs(distPtLine-distPtLine2)>(1e-6)
           warning('Distance of point to line differs by %.1e!', ...
               abs(distPtLine-distPtLine2)); 
        end
    end
end