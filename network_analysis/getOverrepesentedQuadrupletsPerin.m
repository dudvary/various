function [Connectivity,motifID] = getOverrepesentedQuadrupletsPerin()
% Returns the 27 overrepresented Quadruplet motifs from 
% Perin, R., Berger, T. K., & Markram, H. (2011). A synaptic organizing 
% principle for cortical neuronal groups. Proceedings of the National 
% Academy of Sciences of the United States of America, 108(13), 5419–5424. 
% https://doi.org/10.1073/pnas.1016051108
%
% Taken from Fig. S2 "Overexpressed four-neuron connectivity patterns"
% Node IDs:
%    1
%   /|\
% 4 -|- 2
%   \|/
%    3
% Output: 
% - Connectivity [27 x 4 x 4]
%       27 motifs with 4 x 4 connectivity matrix 
% - motifID [1 x 27] code from Figure S2
    
    Connectivity = false(27,4,4); 
    
    % code = 14
    m = false(4,4); 
    m(1,4) = 1; % 1 is connected to 4 
    m(2,1) = 1; 
    m(2,4) = 1; 
    m(2,3) = 1; 
    Connectivity(1,:,:) = m; 

    % code = 15
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,4) = 1; 
    m(1,2) = 1; 
    m(2,1) = 1; 
    Connectivity(2,:,:) = m; 

    % code = 18
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,4) = 1; 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(2,3) = 1; 
    Connectivity(3,:,:) = m; 

    % code = 21
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,4) = 1; 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(2,3) = 1; 
    m(1,3) = 1; 
    Connectivity(4,:,:) = m; 

    % code = 47
    m = false(4,4); 
    m(1,2) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1; 
    m(1,3) = 1; 
    Connectivity(5,:,:) = m; 

    % code = 48
    m = false(4,4); 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1; 
    m(1,3) = 1; 
    Connectivity(6,:,:) = m; 

    % code = 55
    m = false(4,4); 
    m(2,1) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,4) = 1;
    m(2,3) = 1; 
    m(1,3) = 1; 
    Connectivity(7,:,:) = m; 

    % code = 56
    m = false(4,4); 
    m(1,2) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,4) = 1;
    m(1,3) = 1; 
    Connectivity(8,:,:) = m; 

    % code = 59
    m = false(4,4); 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1;
    m(1,3) = 1; 
    Connectivity(9,:,:) = m; 

    % code = 67
    m = false(4,4); 
    m(4,1) = 1; 
    m(2,1) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1;
    m(1,3) = 1; 
    Connectivity(10,:,:) = m; 

    % code = 78
    m = false(4,4); 
    m(4,3) = 1; 
    m(2,1) = 1; 
    m(2,4) = 1; 
    m(2,3) = 1;
    m(1,3) = 1; 
    Connectivity(11,:,:) = m; 

    % code = 82
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,4) = 1; 
    m(1,3) = 1; 
    m(2,3) = 1;
    m(4,3) = 1; 
    Connectivity(12,:,:) = m; 

    % code = 86
    m = false(4,4); 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(2,4) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(4,3) = 1; 
    m(2,3) = 1; 
    Connectivity(13,:,:) = m; 

    % code = 88
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,1) = 1; 
    m(2,4) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(4,3) = 1; 
    m(2,3) = 1; 
    Connectivity(14,:,:) = m; 

    % code = 89
    m = false(4,4); 
    m(1,4) = 1; 
    m(1,2) = 1; 
    m(2,4) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(4,3) = 1; 
    m(2,3) = 1; 
    Connectivity(15,:,:) = m; 

    % code = 106
    m = false(4,4); 
    m(1,4) = 1; 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(4,2) = 1;
    m(2,3) = 1; 
    m(3,2) = 1; 
    Connectivity(16,:,:) = m; 

    % code = 117
    m = false(4,4); 
    m(1,4) = 1; 
    m(1,2) = 1; 
    m(2,4) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    m(2,3) = 1; 
    Connectivity(17,:,:) = m; 

    % code = 148
    m = false(4,4); 
    m(4,1) = 1; 
    m(1,2) = 1; 
    m(2,1) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    Connectivity(18,:,:) = m; 

    % code = 153
    m = false(4,4); 
    m(4,1) = 1; 
    m(1,2) = 1; 
    m(2,4) = 1; 
    m(4,2) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    m(2,3) = 1; 
    Connectivity(19,:,:) = m; 

    % code = 168
    m = false(4,4); 
    m(1,4) = 1; 
    m(2,1) = 1; 
    m(2,4) = 1; 
    m(4,3) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    Connectivity(20,:,:) = m; 

    % code = 171
    m = false(4,4); 
    m(1,4) = 1; 
    m(1,2) = 1; 
    m(2,4) = 1; 
    m(4,3) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    Connectivity(21,:,:) = m; 

    % code = 172
    m = false(4,4); 
    m(2,1) = 1; 
    m(1,2) = 1; 
    m(1,4) = 1; 
    m(2,4) = 1;
    m(1,3) = 1; 
    m(3,2) = 1; 
    m(4,3) = 1; 
    Connectivity(22,:,:) = m; 

    % code = 180
    m = false(4,4); 
    m(2,1) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1;
    m(3,2) = 1; 
    m(1,3) = 1; 
    m(4,3) = 1; 
    Connectivity(23,:,:) = m; 

    % code = 184
    m = false(4,4); 
    m(1,2) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,3) = 1;
    m(3,2) = 1; 
    m(1,3) = 1; 
    m(4,3) = 1; 
    Connectivity(24,:,:) = m; 

    % code = 186
    m = false(4,4); 
    m(1,2) = 1; 
    m(1,4) = 1; 
    m(4,2) = 1; 
    m(2,4) = 1;
    m(3,2) = 1; 
    m(2,3) = 1; 
    m(1,3) = 1; 
    m(4,3) = 1; 
    Connectivity(25,:,:) = m; 

    % code = 188
    m = false(4,4); 
    m(2,1) = 1; 
    m(1,4) = 1; 
    m(4,1) = 1; 
    m(2,4) = 1;
    m(4,3) = 1; 
    m(1,3) = 1; 
    m(3,2) = 1; 
    Connectivity(26,:,:) = m; 

    % code = 189
    m = false(4,4); 
    m(2,1) = 1; 
    m(1,4) = 1; 
    m(4,1) = 1; 
    m(4,3) = 1; 
    m(1,3) = 1; 
    m(3,2) = 1; 
    m(2,3) = 1; 
    Connectivity(27,:,:) = m; 

    % Check whether all are unique
    t = nan(size(Connectivity,1),numel(m));
    for i = 1:size(Connectivity,1)
       m = squeeze(Connectivity(i,:,:));
       t(i,:) = m(:); 
    end

    if size(unique(t,'rows'),1)~=size(Connectivity,1)
       error('Motifs are not all unique!'); 
    end
    
    % code according to Fig S2
    motifID = [14 15 18 21 47 48 55 56 59 67 78 82 86 88 89 ...
                106 117 148 153 168 171 172 180 184 186 188 189]; 
    
end