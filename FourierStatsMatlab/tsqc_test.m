function output = tsqc_test(x,varargin)

% tsqc_test: function to calculate the T-squared circ statistic from Victor & Mast (1991)
% the inputs must be Nx2 or 2xN matrices of numbers
% for a one-sample test, mu is an optional vector to which the data are compared
% for a two-sample test, x and y are the data from the two conditions
% if paired=1, the matrices (x and y) must be the same size
% the output is a structure containing the test results
% this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

y = [];
paired = 0;
mu = [];
if nargin>1
    y = varargin{1};
    if nargin>2
        paired = varargin{2};
        if nargin>3
            mu = varargin{3};
        end
    end
end

d = size(x);
if (d(1)<d(2))
    x = x';
end

if (~isempty(y))
d = size(y);
if (d(1)<d(2))
    y = y';
end
end

if (~isempty(mu))
    for n = 1:length(mu)
        x(:,n) = x(:,n) - mu(n);
    end
end

if (isempty(paired))
    paired = 0;
end

if isreal(x)
    x = complex(x(:,1),x(:,2));
end
if (~isempty(y))
    if isreal(y)
        y = complex(y(:,1),y(:,2));
    end
end

if (isempty(y))
    method = 'One-sample T-squared circ test';
end
if (~isempty(y) && paired==1)
    method = 'Paired samples T-squared circ test';
end
if (~isempty(y) && paired==0)
    method = 'Independent samples T-squared circ test';
end

if (paired==1)
    if (~isempty(y))
        x = x - y;
    end
end

df1 = 2;
if (paired==1 || isempty(y))
    % paired or one-sample version of the test
    
    nobs = length(x);
    cohmean = mean(x);
    absmean = abs(cohmean);
    
    displacement = abs(x - cohmean);
    diffsumsq = sum(displacement.^2);
    tsqc = (nobs-1) .* (absmean.^2) ./ diffsumsq;
    Fratio = nobs*tsqc;
    df2 = nobs*2 - 2;
    pval = fcdf(Fratio,df1,df2,'upper');
    
else
    % independent samples version of the test
    
    nobs1 = length(x);
    cohmean1 = mean(x);
    absmean1 = abs(cohmean1);
    displacement1 = abs(x - cohmean1);
    
    nobs2 = length(y);
    cohmean2 = mean(y);
    absmean2 = abs(cohmean2);
    displacement2 = abs(y - cohmean2);
    
    meandiff = cohmean1 - cohmean2;
    absdiff = abs(meandiff);
    diffsumsq = sum(displacement1.^2) + sum(displacement2.^2);
    
    tsqc = (nobs1 + nobs2 - 2) .* (absdiff.^2) ./ diffsumsq;
    Fratio = ((nobs1*nobs2)./(nobs1 + nobs2)).*tsqc;
    df2 = (2*nobs1 + 2*nobs2 - 4);
    pval = fcdf(Fratio,df1,df2,'upper');
    
end

pval = min(pval,1);

output.tsqc = tsqc;
output.Fratio = Fratio;
output.df1 = df1;
output.df2 = df2;
output.pval = pval;
output.method = method;

end