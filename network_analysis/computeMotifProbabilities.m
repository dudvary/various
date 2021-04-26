function motif = computeMotifProbabilities(p)
%% Compute probability of each motif given a fixed (average) Connection Probability
% Motif order as in
% Egger, R., Dercksen, V. J., Udvary, D., Hege, H., & Oberlaender, M. (2014).
% Generation of dense statistical connectomes from sparse morphological data. 
% Frontiers in Neuroanatomy, 8(November), 1–18. 
% https://doi.org/10.3389/fnana.2014.00129
% 
% 64 motifs in total, reduced to 16
% sum(motif) = 1
% Input:
% - (avg) connection probability
% Output:
% - motif: 1 x 16 with probability of each motif occuring

    if length(p)~=1
       error('only one value as input allowed'); 
    end

    if (p<0 || p>1)
        error('probability has to be between 0 and 1')
    end

    motif = nan(1,16);
    motif(1) = p^6;
    motif(2) = ((1-p)^3*p^3)*2;
    motif(3) = ((1-p)^1*p^5)*3*2;
    motif(4) = ((1-p)^2*p^4)*3*2;
    motif(5) = ((1-p)^2*p^4)*3;
    motif(6) = ((1-p)^2*p^4)*3;
    motif(7) = ((1-p)^3*p^3)*3*2;
    motif(8) = ((1-p)^2*p^4)*3;
    motif(9) = ((1-p)^4*p^2)*2*3;
    motif(10) = ((1-p)^3*p^3)*2*3;
    motif(11) = ((1-p)^3*p^3)*2*3;
    motif(12) = ((1-p)^4*p^2)*3;
    motif(13) = ((1-p)^4*p^2)*3;
    motif(14) = ((1-p)^4*p^2)*3;
    motif(15) = ((1-p)^5*p)*6; 
    motif(16) = (1-p)^6; 

end