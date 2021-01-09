function output = CI_test(data,alpha)
% CI.test: function to calculate the condition index of bivariate data
% Inputs:
%   data is an N x 2 matrix of bivariate (x,y) observations, or a vector of N complex values
%   alpha is the criterion for statistical significance, set to a default of 0.05
% Output:
%   a structure with fields: CI (the condition index), N (sample size),
%   criticalCI (the critical index, given alpha), pval (the p-value)
% the condition index is square root of the the ratio of eigenvalues (eigenvector lengths), calculated longest/shortest
% the index is then compared for significance with an expected distribution function
% significant tests (p < alpha) violate the assumptions of the T-squared-circ and ANOVA-squared-circ tests
% see Baker (2021) for further details
% this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

if (isempty(alpha))
    alpha = 0.05;
end

d = size(data);
if (d(1)<d(2))
    data = data';
end

if (~isreal(data))
    data = [real(data) imag(data)];
end

C = cov(data);
N = length(data);
eigVal = eig(C);
CI = sqrt(max(eigVal)/min(eigVal));

cilist = 1:0.001:100;

pdffunction = ((N-2).*(2.^(N-2))) .* ((cilist.^2 - 1)./((cilist.^2+1).^(N-1))) .* (cilist.^(N-3));
cdfinverse = 1 - (cumsum(pdffunction)./sum(pdffunction));
criticalCI = min(cilist(find(cdfinverse<alpha)));

pval = 0;
indices = find(cilist>=CI);
if (~isempty(indices))
    pval = cdfinverse(indices(1));
end

output.CI = CI;
output.N = N;
output.criticalCI = criticalCI;
output.pval = pval;
    
end