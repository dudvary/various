%% Perform a Kalmagorov-Smirnov Test based on Histograms (not samples as implemented in kstest2)
% Checks whether there's a significant difference between two histograms
% (two-sided test)
% Input: 
% - hist1: histogram of samples1
% - hist2: histogram of samples2
% - alpha: significant level (0.05)
% Output:
% - H: [0 1] for rejection H0 (same Distribution)
% - pValue: p-value
% - KSstatistic: Max. Distance between empirical CDFs
function [H, pValue, KSstatistic] = kstest2_hist(hist1,hist2,alpha)

    % H0: Same Distribution
    % if h = 1, reject H0
    n1 = sum(hist1);                        % Sample Size 1
    n2 = sum(hist2);                        % Sample Size 2

    s1 = cumsum(hist1)./n1;                 % Empirical CDF 1
    s2 = cumsum(hist2)./n2;                 % Empirical CDF 2

    KSstatistic = max(abs(s1-s2));          % Max. Difference of Histograms

    % Copied from kstest2
    n = n1 * n2 /(n1 + n2);
    lambda = max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);
    j =  (1:101)';
    pValue = 2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
    pValue = min(max(pValue, 0), 1);

    H  =  (alpha >= pValue);
end