function output = pairwisemahal(data,grouping)

% calculates the pairwise (group to group) Mahalanobis distance between group means
%
% inputs:
% data - can be either a vector of N complex numbers, a matrix of Nxm complex numbers, or an Nx2 matrix of real numbers
%        if a vector of complex numbers or matrix of real numbers, the 'grouping' variable is also required
%        if a matrix of complex numbers is supplied, the columns are interpreted as the groups
% grouping - a variable containing group membership information for the data
%
% output structure (each variable is a subfield):
% D - an mxm matrix of Mahalanobis distances between the group means
%     note that the Mahalanobis distance (D) is returned, not D^2
% groups - a list of group names in the order corresponding to the rows and
%          columns of the matrix D
%
% this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

d = size(data);
if (d(1)<d(2))
    data = data';
end
d = size(data);

% convert data to an Nxm matrix of real numbers
if (~isreal(data))
    if (min(d)==1)  % it is a vector
        realdata = [real(data) imag(data)];         % convert to real
    else            % it is a matrix
        ngroups = d(2);
        realdata = zeros(d(1),2);
        for n = 1:ngroups
            realdata(((n-1)*d(1)+1):(n*d(1)),:) = [real(data(:,n)) imag(data(:,n))];
            grouping(((n-1)*d(1)+1):(n*d(1))) = n;
        end
        
    end
else
    realdata = data;
end

grouplabels = unique(grouping);
m = length(grouplabels);

for n = 1:m
    centroids(n,:) = mean(realdata(find(grouping==grouplabels(n)),:));
end

for n1 = 1:m
    for n2 = 1:m
        CA = cov(realdata(find(grouping==grouplabels(n1)),:)); % covariance matrix for group A
        CB = cov(realdata(find(grouping==grouplabels(n2)),:)); % covariance matrix for group B
        N1 = length(find(grouping==grouplabels(n1)));          % sample size for group A
        N2 = length(find(grouping==grouplabels(n2)));          % sample size for group B
        pwcov = ((N1-1)*CA + (N2-1)*CB)/(N1+N2-2);            % pooled covariance matrix for the two groups
        meandiff = centroids(n1,:) - centroids(n2,:);         % difference in means
        distances(n1,n2) = sqrt(meandiff * inv(pwcov) * meandiff'); % Mahalanobis distance
    end
end

output.D = distances;
output.groups = grouplabels;

end