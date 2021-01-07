function output = tsqh_test(x,y,paired,mu)

% tsqh.test: function to calculate a different variants of the T-squared test of Hotelling (1931)
% the inputs must be Nx2 or 2xN matrices of numbers
% for a one-sample test, mu is an optional vector to which the data are compared
% for a two-sample test, x and y are the data from the two conditions
% if paired=TRUE, the matrices (x and y) must be the same size

d = size(x);
if (d(1)<d(2))
    x = x';
end
d = size(y);
if (d(1)<d(2))
    y = y';
end

if (~isempty(mu))
    for n = 1:length(mu)
        x(:,n) = x(:,n) - mu(n);
    end
end

if ~isreal(x)
    x = [real(x); imag(x)]';
end
if (~isempty(y))
    if ~isreal(y)
        y = [real(y); imag(y)]';
    end
end

if (isempty(y))
    method = 'One-sample T-squared test';
end
if (~isempty(y) && paired==1)
    method = 'Paired samples T-squared test';
end
if (~isempty(y) && paired==0)
    method = 'Independent samples T-squared test';
end

if (paired==1)
    if (~isempty(y))
        x = x - y;
    end
end


if (paired==1 || isempty(y))
    % paired or one-sample version of the test

    xbar = mean(x);
    s = size(x);
    m = s(2);
    N = s(1);
    C = cov(x);
    Cinv = inv(C);
    tsq = N * xbar * Cinv * xbar';
    Fratio = tsq .* (N - m)./(m.*(N-1));
    df1 = m;
    df2 = N-m;
    pval = fcdf(Fratio,df1,df2,'upper');
 
else
    % independent samples version of the test
    
    xbarA = mean(x);
    s = size(x);
    m = s(2);
    N1 = s(1);
    
    xbarB = mean(y);
    s = size(y);
    N2 = s(1); 
    
    CA = cov(x);
    CB = cov(y);
    
    meandiff = xbarA - xbarB;
    
    C = ((N1-1)*CA + (N2-1)*CB)/(N1+N2-m);
    Cinv = inv(C);
    
    tsq = ((N1*N2)/(N1+N2)) * meandiff * Cinv * meandiff';
    Fratio = tsq * (N1+N2-m-1)/(m*(N1+N2-2));
    df1 = m;
    df2 = N1+N2-m-1;
    
    pval = fcdf(Fratio,df1,df2,'upper');
    
end

pval = min(pval,1);

output.tsqc = tsq;
output.Fratio = Fratio;
output.df1 = df1;
output.df2 = df2;
output.pval = pval;
output.method = method;

end