%% Calculates on which side of a line (defined by v1 and v2) a point pt
% lies
% Input:
% - pt: point of interest 
% - v1: vector1 that spans the line
% - v2: vector2 that spans the line
% Output:
% - s: either -1,0,or 1
function s = sideOfLine(pt,v1,v2)

    if (length(v1)==2 && length(v2)==2 && length(pt)==2)
        s = (v2(1)-v1(1))*(pt(2)-v1(2))-(v2(2)-v1(2))*(pt(1)-v1(1));
        s = sign(s); 
    else
        error(['Input has to be 2-element arrays!']);
    end
end