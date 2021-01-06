function output = anovacirc_test(data,group,participant)

% anovacirc.test: two-dimensional analysis of variance using complex data
% an extension of the logic of the T-squared circ statistic of Victor & Mast (1991)
% this is a one-way implementation of the test - if participant IDs are supplied a repeated measures version is conducted
%
% Inputs--
% data:  this can either be a vector of complex numbers, or a matrix
%        if it is a matrix, either the first two columns are the x and y (real and imaginary) values of the DV
%        if there are further columns, these are treated as the group labels and participant IDs
% group: condition labels indicating the level of the independent variable that each data point corresponds to
%        this is an optional input, but if it is supplied it supercedes values from the matrix
% participant: variable storing participant/subject IDs (i.e. the random factor) for repeated measures analysis
%        this is an optional input, but if it is supplied it supercedes values from the matrix
%
% see Baker (2021) for further details

grouplabels = [];
participantlabels = [];

s = size(data);

if (~isreal(data(1)))
    datavals = data;
else
    datavals = complex(data(:,1),data(:,2));
    
    if (s(2)>2)
        grouplabels = data(:,3);
    end
    if (s(2)>3)
        participantlabels = data(:,4);
    end
end

if (~isempty(group))
    grouplabels = group;
end
if (~isempty(participant))
    participantlabels = participant;
end

factorlist = unique(grouplabels);

grandmean = mean(datavals);

if (isempty(participantlabels))
    % if no participant labels have been supplied run a between-subject ANOVA-circ test
    
    SSM = 0;
    allresiduals = zeros(1,length(datavals));
    for n = 1:length(factorlist)
        groupmeans(n) = mean(datavals(find(grouplabels==factorlist(n))));
        SSM = SSM + length(find(grouplabels==factorlist(n))) .* abs(groupmeans(n)-grandmean).^2;
        allresiduals(find(grouplabels==factorlist(n))) = abs(datavals(find(grouplabels==factorlist(n))) - groupmeans(n));
    end
    
    dfM = 2*(length(factorlist)-1);
    SSR = sum(allresiduals.^2);
    dfR = 2*(length(datavals)-length(factorlist));
    MSM = SSM/dfM;
    MSR = SSR/dfR;
    Fratio = MSM/MSR;
    pval = fcdf(Fratio,dfM,dfR,'upper');
    pval = min(pval,1);
    
    output.Fratio = Fratio;
    output.pval = pval;
    output.SSM = SSM;
    output.SSR = SSR;
    output.dfM = dfM;
    output.dfR = dfR;
    output.MSM = MSM;
    output.MSR = MSR;
    output.method = 'Between-subjects ANOVA2-circ';
    
else
    % participant labels have been supplied, run a repeated measures ANOVA-circ test
    
    participantlist = unique(participantlabels);
    
    SSW = 0;
    for n = 1:length(participantlist)
        SSW = SSW + sum(abs(datavals(find(participantlabels==participantlist(n))) - mean(datavals(find(participantlabels==participantlist(n))))).^2);
    end
    dfW = 2*(length(participantlist)*(length(factorlist)-1));
    
        SSM = 0;
    for n = 1:length(factorlist)
        groupmeans(n) = mean(datavals(find(grouplabels==factorlist(n))));
        SSM = SSM + length(find(grouplabels==factorlist(n))) .* abs(groupmeans(n)-grandmean).^2;
    end
    dfM = 2*(length(factorlist)-1);
    
    SSR = SSW - SSM;
    dfR = dfW - dfM;
    MSM = SSM/dfM;
    MSR = SSR/dfR;
    Fratio = MSM/MSR;
    pval = fcdf(Fratio,dfM,dfR,'upper');
    pval = min(pval,1);
    
    output.Fratio = Fratio;
    output.pval = pval;
    output.SSW = SSW;
    output.SSM = SSM;
    output.SSR = SSR;
    output.dfW = dfW;
    output.dfM = dfM;
    output.dfR = dfR;
    output.MSM = MSM;
    output.MSR = MSR;
    output.method = 'Repeated measures ANOVA2-circ';
   
end

end