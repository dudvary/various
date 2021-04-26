function [R,R4x4] = computeTransformationMatrix(v,u)
% compute Transformation Matrix (3x3)
% rotates v to u so that
% v * R = u
% Transformation Matrix (4x4)
% [v 1] * R4x4 = u4x4; 
% Input:
% - v: vector to be rotated (1x3)
% - u: final vector to be rotated to (1x3)
%       (optional) default: 0 0 1
% Output:
% - R: transformationMatrix 3x3 so that v * R = u
% - R4x4: transformationMatrix 4x4 so that [v 1] * R4x4 = u4x4
%       for scaling and translation

    if nargin==1
        u = [0 0 1];
    end
    
    if length(v)~=3 || length(u)~=3
       error('vector has to be 1x3'); 
    end

    v = v./norm(v);
    u = u./norm(u); 
    w = cross(u,v)/norm(cross(u,v)); 
    U = [u' w' cross(u,w)']; 
    V = [v' w' cross(v,w)']; 
    R = V*U'; 
    
    if nargout>1
        R4x4 = eye(4,4);
        R4x4(1:3,1:3) = R;
    end

end

