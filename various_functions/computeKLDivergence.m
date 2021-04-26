function [KLdiv] = computeKLDivergence(P,Q)
%% Calculate Kullback-Leibler-Divergence KL(P||Q)
% between two discrete Probability Distributions
% Input
% - P: probability function [observed/true distribution]
% - Q: probability function [function]
% Output:
% - KLdiv: Kullback-Leibler-Divergence
%       -> information lost, the higher the more information is lost
%       -> KLdiv = 0, both are identical
%       -> expected number of bits required when using a code based on Q 
%       rather than P -> expected number of extra bits that must be 
%       transmitted to identify a value x drawn from X, if a code is used 
%       corresponding to probability distribution Q, rather than the 
%       "true" distribution P.

    if length(P) ~= length(Q)
        error('P and Q have to be same lenght!');
    end
    
    idx = find(Q>1e-10 & P>1e-10); 
    
    if length(idx)<2
       warning('Less than 2 valid bins found. KL set to NaN');
       KLdiv = nan;
       return; 
    end
    
    Q = Q(idx);
    P = P(idx); 

    if sum(Q)~=1
        Q = Q./sum(Q);
    end
    
    if sum(P)~=1
        P = P./sum(P);
    end
    
    KLdiv = sum(P.*log(P./Q)); 
    
end
