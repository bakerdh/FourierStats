function output = analysecplx(data, group, participant)

% analysecplx: wrapper function that automatically decides which statistical test to run for a given data set
%
% Inputs--
% data:  this can either be a vector of complex numbers, or a matrix
%        if it is a matrix, the first two columns are the x and y (real and imaginary) values of the DV
%        if there are further columns, these are treated as the group labels and participant IDs
% group: condition labels indicating the level of the independent variable that each data point corresponds to
%        this is an optional input, but if it is supplied it supercedes values from the matrix
% participant: variable storing participant/subject IDs (i.e. the random factor) for repeated measures analysis
%        this is an optional input, but if it is supplied it supercedes values from the matrix
%
% the function first runs the condition index test for each level of the IV
% if the condition index test is non-significant for all levels, a T-squared-circ test is run for one- and two-sample designs, or an ANOVA-squared-circ test otherwise
% if the condition index test is significant for any level, a T-squared test or MANOVA is run instead
% if participant labels are included, a repeated measures test is conducted
% the Mahalanobis distance effect size statistic (D) is also calculated

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

if (isempty(grouplabels))
    grouplabels = ones(size(complexdata));
end
isRM = 1;
if (isempty(participantlabels))
    isRM = 0;
end

factorlist = unique(grouplabels);
ngroups = length(factorlist);
nobservations = length(complexdata);

disp(strcat('Design has ',num2str(ngroups),' levels, with ',numstr(nobservations/ngroups),' observations per level'));
if (isRM==1)
    disp('Design is repeated measures');
end


% run the condition index test for each group
CIresults = zeros(ngroups,4);
for n = 1:ngroups
    i = find(grouplabels==factorlist(n));
    temp = CI_test(complexdata(i));
    CIresults = [temp.CI temp.N temp.criticalCI temp.pval];
end
nsigCI = length(find(CIresults(:,4)<(0.05/ngroups)));
disp(strcat('The condition index was significant for ',num2str(nsigCI),' out of ',num2str(ngroups),' levels'));
assumptionsmet = 1;
if (nsigCI>0)
    assumptionsmet = 0;
end
output.ngroups = ngroups;
output.nsigCI = nsigCI;


% run T-squared-circ
if (assumptionsmet==1 && ngroups<3)
    if (ngroups==1)
        disp('Running a one-sample T-squared-circ test');
        output.testtype = 'One-sample T-squared-circ test';
        results = tsqc_test(complexdata);
    else
        if (isRM==0)
            disp('Running an independent T-squared-circ-test');
            output.testtype = 'Independent T-squared-circ test';
        else
            disp('Running a repeated measures T-squared-circ-test');
            output.testtype = 'Repeated measures T-squared-circ test';
        end
        dataA = complexdata(find(grouplabels==factorlist(1)));
        dataB = complexdata(find(grouplabels==factorlist(2)));
        results = tsqc_test(dataA,dataB,isRM);
        disp(strcat('The test statistic T^2_circ = ',num2str(results.tsqc,2)));
        disp(strcat('The equivalent F-ratio with ',num2str(results.df1),' and ',num2str(results.df2),' degrees of freedom is F = ',num2str(results.Fratio,2)));
        output.teststat = results.tsqc;
        output.Fratio = results.Fratio;
        output.df1 = results.df1;
        output.df2 = results.df2;
        output.pval = results.pval;
    end
end


% run T-squared
if (assumptionsmet==0 && ngroups<3)
    if (ngroups==1)
        disp('Running a one-sample T-squared test');
        output.testtype = 'One-sample T-squared test';
        results = tsq1.test(complexdata);
        disp(strcat('The test statistic T^2 = ',num2str(results.tsq,2)));
        disp(strcat('The equivalent F-ratio with ',num2str(results.df1),' and ',num2str(results.df2),' degrees of freedom is F = ',num2str(results.Fratio,2)));
        output.teststat = results.tsq;
        output.Fratio = results.Fratio;
        output.df1 = results.df1;
        output.df2 = results.df2;
        output.pval = results.pval;
    end
    if (ngroups>1)
        dataA = complexdata(find(grouplabels==factorlist(1)));
        dataB = complexdata(find(grouplabels==factorlist(2)));
        if (isRM==0)
            disp('Running an independent T-squared-test');
            output.testtype = 'Independent T-squared test';
            
            
            % need to run an independent T-squared test here
            
            
        else
            disp('Running a repeated measures T-squared-test');
            output.testtype = 'Repeated measures T-squared test';
            results = tsq1.test(dataA-dataB);
            disp(strcat('The test statistic T^2 = ',num2str(results.tsq,2)));
            disp(strcat('The equivalent F-ratio with ',num2str(results.df1),' and ',num2str(results.df2),' degrees of freedom is F = ',num2str(results.Fratio,2)));
            output.teststat = results.tsq;
            output.Fratio = results.Fratio;
            output.df1 = results.df1;
            output.df2 = results.df2;
            output.pval = results.pval;
        end
    end
end


% run ANOVA-squared-circ
if (assumptionsmet==1 && ngroups>2)
    if (isRM==0)
        disp('Running an independent ANOVA-squared-circ test')
        output.testtype = 'Independent ANOVA-squared-circ test';
    else
        disp('Running a repeated measures ANOVA-squared-circ test');
        output.testtype = 'Repeated measures ANOVA-squared-circ test';
    end
    results = anovacirc_test(complexdata,grouplabels,participantlabels);
    disp(strcat('The F-ratio with ',num2str(results.dfM),' and ',num2str(results.dfR),' degrees of freedom is F = ',num2str(results.Fratio,2)));
    output.Fratio = results.Fratio;
    output.df1 = results.dfM;
    output.df2 = results.dfR;
    output.pval = results.pval;
end


% run MANOVA
if (assumptionsmet==0 && ngroups>2)
if (isRM==0)
    disp('Running an independent MANOVA');
    output.testtype = 'Independent MANOVA (test statistic is Wilks Lambda)';
    data = [real(complexdata); imag(complexdata)]';
    [D,P,STATS] = manova1(data,grouplabels);
    output.teststat = STATS.lambda;
    output.pval = P(1);
else
    % repeated measures MANOVA not properly implemented
    disp('Data set requires repeated measures MANOVA, which is not available in Matlab');
    disp('Consider using the R version of this toolbox (or SPSS) for this analysis');
end
    
end


if (output.pval<0.05)
    disp('The test was significant at p < 0.05');
else
    disp('The test was not significant at p < 0.05');
end


temp = [real(complexdata); imag(complexdata)]';
if (ngroups==1)
    % calculate pointwise Mahalanobis distance relative to the origin
    D = sqrt(mahal([0,0],temp));
else
    % calculate pairwise Mahalanobis distance between each pair of conditions
    % this actually uses the manova1 function
    [D,P,STATS] = manova1(temp,grouplabels);   % the gmdist field contains pairwise Mahalanobis statistics
    D = max(STATS.gmdist);   % return the largest effect size
end
disp(strcat('The effect size (Mahalanobis distance) is D = ',num2str(D,2)));
output.D = D;


end