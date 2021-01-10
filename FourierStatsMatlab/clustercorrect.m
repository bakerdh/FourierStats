function clustout = clustercorrect(datax,varargin)

% clustercorrect: implements non-parametric cluster correction method described by Maris & Oostenveld (2007) J Neurosci Methods, 164: 177-190, doi: 10.1016/j.jneumeth.2007.03.024
% this version can deal with complex data using the T-squared or T-squared-circ statistics
%
% function call: clustercorrect(datax, datay, adjacencymatrix, testtype, paired, clustformthresh, clustthresh, nresamples)
% The input datax should be an N (participants) x m (sensors/timepoints/locations) matrix
% The optional input datay can either be a matrix of the same dimensions, or N2 x m (for unpaired designs with unbalanced samples), or a single value to compare to datax
%
% Permitted values of the testtype input are 1: t-test, 2: Hotelling's T-squared, 3: Victor & Mast's T-squared-circ
% Complex values are required for options 2 & 3. If complex values are passed in for test 1 (t-test), they will be converted to amplitudes
%
% if the measures are from adjacent time points, an adjacency matrix will be automatically constructed
% if the measures are from different spatial locations (e.g. an electrode montage) you should calculate an adjacency matrix yourself
% this should be an m x m matrix (where m is the number of sensors), where values are zero for non-adjacent sensor pairs, and 1 for adjacent pairs
%
% outputs an object containing the indices of all significant clusters
% this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

datay = [];
adjacencymatrix = [];
testtype = 3;
paired = 0;
clustformthresh = 0.05;
clustthresh = 0.05;
nresamples = 1000;

if nargin>1
    datay = varargin{1};
    if nargin>2
        adjacencymatrix = varargin{2};
        if nargin>3
            testtype = varargin{3};
            if nargin>4
                paired = varargin{4};
                if nargin>5
                    clustformthresh = varargin{5};
                    if nargin>6
                        clustthresh = varargin{6};
                        if nargin>7
                            nresamples = varargin{7};
                        end
                    end
                end
            end
        end
    end
end


clustout = [];
sigclustercount = 0;
clusterpout = [];
clusterpoints = {};
s = size(datax);
N = s(1);  % number of observations (e.g. participants)
m = s(2);  % number of measurements (e.g. sensors, voxels, timepoints)

% if complex data have been passed in and t-test selected, convert to amplitudes
if (~isreal(datax) && testtype==1)
    datax = abs(datax);
end
if (~isreal(datay) && testtype==1)
    datay = abs(datay);
end

% if two data sets have been passed in, note the size of the second matrix
if (~isempty(datay))
    N2 = length(datay);
end

% if it's a paired design, subtract the two data sets and delete the second one
if (paired==1 && ~isempty(datay))
    datax = datax - datay;
    datay = [];
end

if (isempty(adjacencymatrix))
    % if no adjacency matrix has been passed in, we assume that subsequent observations are adjacent
    adjacencymatrix = zeros(m,m);
    for n = 1:(m-1)
        adjacencymatrix(n,n+1) = 1;
    end
    for n = 2:m
        adjacencymatrix(n,n-1) = 1;
    end
end

% loop through all elements (sensors/time points) and calculate the test statistic and p-value

for n = 1:m
    
    if (isempty(datay))
        % one-sample or paired test
        
        switch testtype
            case 1  % univariate t-test
                [h,p,ci,output] = ttest(datax(:,n));
                allp(n) = p;
                allt(n) = output.tstat;
            case 2  % one-sample Hotelling's t-squared
                output = tsqh_test(datax(:,n));
                allp(n) = output.pval;
                allt(n) = output.tsq;
            case 3  % one-sample T-squared-circ
                output = tsqc_test(datax(:,n));
                allp(n) = output.pval;
                allt(n) = output.tsqc;
        end
        
    else
        % independent samples test
        
        switch testtype
            case 1  % independent univariate t-test
                [h,p,ci,output] = ttest2(datax(:,n),datay(:,n));
                allp(n) = p;
                allt(n) = output.tstat;
            case 2  % two sample Hotelling's t-squared
                output = tsqh_test(datax(:,n),datay(:,n));
                allp(n) = output.pval;
                allt(n) = output.tsq;
            case 3  % two-sample T-squared-circ
                output = tsqc_test(datax(:,n),datay(:,n));
                allp(n) = output.pval;
                allt(n) = output.tsqc;
        end
        
    end
    
end
allp(find(isnan(allp))) = 1;

% now generate a list of clusters of adjacent significant elements

hvect = zeros(1,m);
hvect(allp<clustformthresh) = 1;
hmat = repmat(hvect,[m 1]);
hmat2 = hmat .* hmat';
clustprod = hmat2 .* adjacencymatrix;

clusterlist = {};
nclusts = 0;
for n = 1:m
    temp = sum(clustprod(n,:));
    if temp>0
        i = find(clustprod(n,:)>0);
        tempclust = [n i];
        clustprod(n,:) = 0;
        clustprod(:,n) = 0;
        cnt = 1;
        while cnt<length(tempclust)
            cnt = cnt + 1;
            temp = sum(clustprod(tempclust(cnt),:));
            if temp>0
                i = find(clustprod(tempclust(cnt),:)>0);
                tempclust = [tempclust i];
                clustprod(tempclust(cnt),:) = 0;
                clustprod(:,tempclust(cnt)) = 0;
            end
        end
        nclusts = nclusts + 1;
        clusterlist{nclusts} = unique(tempclust);
    end
end

if nclusts>0
for cc = 1:nclusts
    Cindices = clusterlist{cc};
    sumtvals(cc) = sum(allt(Cindices));
end

i = find(sumtvals==max(sumtvals)); % find the largest cluster
maxcluster = clusterlist{i}; % store the largest cluster

if (length(maxcluster)>1)
    % build a null distribution by permuting signs/group labels
    for n = 1:nresamples
        tsum = 0;
        if (isempty(datay))
            randsigns = (round(rand(N,1))*2)-1;
        else
            randgroups = randperm(N+N2);
        end
        for elcounter = 1:length(maxcluster)
            
            if (isempty(datay))  % one sample or paired test
                tempdataA = datax(:,maxcluster(elcounter)).*randsigns;
                
                switch testtype
                    case 1  % univariate t-test
                        [h,p,ci,output] = ttest(tempdataA);
                        tsum = tsum + output.tstat;
                    case 2  % one-sample Hotelling's t-squared
                        output = tsqh_test(tempdataA);
                        tsum = tsum + output.tsq;
                    case 3  % one-sample T-squared-circ
                        output = tsqc_test(tempdataA);
                        tsum = tsum + output.tsqc;
                end
                
                
            else
                tempdata = [datax(:,maxcluster(elcounter)),datay(:,maxcluster(elcounter))];
                tempdataA = tempdata(randgroups(1:N));
                tempdataB = tempdata(randgroups((N+1):(N+N2)));
                switch testtype
                    case 1  % independent univariate t-test
                        [h,p,ci,output] = ttest2(tempdataA,tempdataB);
                        tsum = tsum + output.tstat;
                    case 2  % two sample Hotelling's t-squared
                        output = tsqh_test(tempdataA,tempdataB);
                        tsum = tsum + output.tsq;
                    case 3  % two-sample T-squared-circ
                        output = tsqc_test(tempdataA,tempdataB);
                        tsum = tsum + output.tsqc;
                end
            end
            nulldist(n) = tsum;
        end
    end
    
    % compare each cluster to the null distribution, retain the significant ones
    
    clusterps = zeros(1,nclusts);
    if (length(nulldist)==nresamples)
        for cc = 1:nclusts
            i = find(abs(nulldist)>abs(sumtvals(cc)));
            if (~isempty(i))
                clusterps(cc) = length(i)/nresamples;
            end
            if (clusterps(cc)<clustthresh)
                sigclustercount = sigclustercount + 1;
                clusterpoints{sigclustercount} = clusterlist{cc};
                clusterpout(sigclustercount) = clusterps(cc);
            end
        end
    end
end
end


clustout.clusterpoints = clusterpoints;
clustout.nclusters = sigclustercount;
clustout.pvals = clusterpout;

end