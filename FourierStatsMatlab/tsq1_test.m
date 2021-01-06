function output = tsq1_test(x,mu)

% tsq1_test: function to calculate a one-sample Hotelling's T-squared statistic 
% the input (x) must be Nx2 or 2xN matrices of numbers
% or can be stored as a single vector of complex numbers
% mu is an optional vector to which the data are compared
% the output is a structure containing the test results

d = size(x);
if (d(1)<d(2))
    x = x';
end

if (~isempty(mu))
    for n = 1:length(mu)
        x(:,n) = x(:,n) - mu(n);
    end
end

if ~isreal(x)
    x = [real(x); imag(x)]';
end

    method = 'One-sample Hotellings t-squared test';

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
    pval = min(pval,1);
    
output.tsq = tsq;
output.Fratio = Fratio;
output.df1 = df1;
output.df2 = df2;
output.pval = pval;
output.method = method;

end