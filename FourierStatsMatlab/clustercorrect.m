function clustout = clustercorrect(datax,datay,adjacencymatrix,testtype,paired,clustformthresh,clustthresh,nresamples)

% clustercorrect: implements non-parametric cluster correction method described by Maris & Oostenveld (2007) J Neurosci Methods, 164: 177-190, doi: 10.1016/j.jneumeth.2007.03.024
% this version can deal with complex data using the T-squared or T-squared-circ statistics
%
% The input datax should be an N (participants) x m (sensors/timepoints/locations) matrix
% The optional input datay can either be a matrix of the same dimensions, or N2 x m (for unpaired designs with unbalanced samples), or a single value to compare to datax
%
% Legal values of the testtype input are 1: t-test, 2: Hotelling's T-squared, 3: Victor & Mast's T-squared-circ
% Complex values are required for options 2 & 3. If complex values are passed in for test 1 (t-test), they will be converted to amplitudes
%
% if the measures are from adjacent time points, an adjacency matrix will be automatically constructed
% if the measures are from different spatial locations (e.g. an electrode montage) you should calculate an adjacency matrix yourself
% this should be an m x m matrix (where m is the number of sensors), where values are zero for non-adjacent sensor pairs, and 1 for adjacent pairs
%
% outputs an object containing the indices of all significant clusters

if (isempty(paired))
    paired = 0;
end
if (isempty(clustformthresh))
    clustformthresh = 0.05;
end
if (isempty(clustthresh))
    clustthresh = 0.05;
end
if (isempty(nresamples))
    nresamples = 1000;
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
        adjacencymatrix(n-1,n) = 1;
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
                output = tsq1_test(datax(:,n));
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
                
            case 3  % two-sample T-squared-circ
                output = tsqc_test(datax(:,n),datay(:,n));
                allp(n) = output.pval;
                allt(n) = output.tsqc;
        end
        
    end
    
end
allp(find(isnan(allp))) = 1;

% now generate a list of clusters of adjacent significant elements

allclusters = {};
clusterstarts = [];
clusterends = [];
nclusters = 0;
incluster = 0;
for n = 1:m
    if (allp(n)<clustformthresh)
        nclusters = nclusters + 1;
        clusterlist = n;
        for n2 = 1:m
            if (adjacencymatrix(n,n2)==1 && allp(n2)<clustformthresh)
                clusterlist = [clusterlist n2];
            end
        end
        allclusters{nclusters} = clusterlist;
    end
end

if (nclusters>0)
    allclusters{nclusters+1} = 0; % fiddle to prevent later errors
    
    % condense the clusters by pooling any with overlapping elements
    
    ccount = 0;
    clustersizes = [];
    condensedclusters = {};
    for n = 1:nclusters
        targetcluster = allclusters{n};
        if (~isempty(targetcluster))
            if (sum(targetcluster)>0)
                for n2 = (n+1):nclusters
                    compcluster = allclusters{n2};
                    if (~isempty(compcluster))
                        if (sum(compcluster)>0)
                            lia = ismember(targetcluster,compcluster);
                            if (sum(lia)>0)
                                targetcluster = [targetcluster, compcluster);
                                allclusters{n2} = 0;
                            end
                        end
                    end
                end
                ccount = ccount + 1;
                condensedclusters{ccount} = unique(targetcluster);
                clustersizes(ccont) = length(unique(targetcluster));
            end
        end
    end
    
    for cc = 1:ccount
        Cindices = condensedclusters{cc};
        sumtvals(cc) = sum(allt(Cindices));
    end
    
    i = find(sumtvals==max(sumtvals)); % find the largest cluster
    maxcluster = condensedclusters{i}; % store the largest cluster
    
    if (length(maxcluster)>1)
        % build a null distribution by permuting signs/group labels
        for n = 1:nresamples
            tsum = 0;
            if (isempty(datay))
                randsigns = (round(rand(1,N))*2)-1;
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
                            output = tsq1_test(tempdataA);
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
                            
                        case 3  % two-sample T-squared-circ
                            output = tsqc_test(tempdataA,tempdataB);
                            tsum = tsum + output.tsqc;
                    end
                end
                nulldist(n) = tsum;
            end
        end
        
        % compare each cluster to the null distribution, retain the significant ones
        
        clusterps = zeros(1,ccount);
        if (length(nulldist)==nresamples)
            for cc = 1:ccount
                i = find(abs(nulldist)>abs(sumtvals(cc)));
                if (~isempty(i))
                    clusterps(cc) = length(i)/nresamples;
                end
                if (clusterps(cc)<clusterthresh)
                    sigclustercount = sigclustercount + 1;
                    clusterpoints{sigclustercount} = condensedclusters{cc};
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