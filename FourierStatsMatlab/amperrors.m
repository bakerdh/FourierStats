function output = amperrors(input,method,quantiles,nresamples)

% amperrors: calculate error bars on the amplitude component of coherently averaged data
% the input is a vector of complex numbers, or an Nx2 matrix of real and imaginary components
% the optional 'method' flag allows the following values:
%  - 'boot': (default) bootstraps the confidence intervals
%  - 'circ': calculates confidence intervals based on a circular bounding region (used when the condition index test is non-significant)
%  - 'ellipse': calculates confidence intervals based on an elliptical bounding region (based on Pei et al. 2017, doi: 10.1016/j.visres.2016.03.010)
%  - 'abs': calculates the standard error (and the mean) using the absolute amplitudes
%  the optional 'quantiles' flag allows the quantile to be set for the confidence intervals
%  the optional 'nresamples' variable allows the number of iterations to be determined for the bootstrapping option (ignored for other methods)
%  typical values are 95 (for 95% confidence intervals), and 68 (for standard errors)
%  the function returns the mean amplitude and the upper and lower error bars

if isempty(quantiles)
    quantiles = 95;
end
if isempty(nresamples)
    nresamples = 10000;
end
if isempty(method)
    method = 'b';
end

s = size(input);
compdata = input;
if (isreal(input))
    compdata = complex(input(:,1),input(:,2));
end
xydata = [real(compdata'); imag(compdata')]';

output.meanamp = abs(mean(compdata));

if (method(1)=='b')
    % bootstrap the error bars on the amplitude using 10000 resamples
    for n = 1:nresamples
        bspop(n) = abs(mean(compdata(randi(length(compdata),1,length(compdata)))));
    end
    bspop = sort(bspop);
    output.lowerCI = bspop(round(nresamples*(0.5*(100-quantiles)/100)));
    output.upperCI = bspop(round(nresamples*(1-(0.5*(100-quantiles)/100))));
end

if (method(1)=='c')
    % pools the variance in the x (real) and y (imaginary) directions
    % calculates the standard error based on this estimate
    sdcirc = sqrt((std(xydata(:,1)).^2 + std(xydata(:,2)).^2)./2);
    se = sdcirc./sqrt(length(compdata));
    if (quantiles==95)
        ci = 1.96*se;
    end
    if (quantiles==68)
        ci = se;
    end
    output.lowerCI = max(output.meanamp - ci,0);  % lower bound at 0
    output.upperCI = output.meanamp + ci;
end


if (method(1)=='e')
    % implements the method of Pei et al. (2017)
    % calculates the bounding ellipse for the data points
    % then uses this to find the nearest and farthest points from the origin
    % these points are used as the lower and upper error bars
    
    angles = linspace(0,2*pi,200);  % angles for ellipse
    ctr = mean(xydata);
    [eigVec,eigenval] = eig(cov(xydata));
    eigVal = 1.96*eigenval/sqrt(length(xydata));
    eigScl = eigVec * sqrt(eigVal);
    xMat = [ctr(1) + eigScl(1,:); ctr(1) - eigScl(1,:)];
    yMat = [ctr(2) + eigScl(2,:); ctr(2) - eigScl(2,:)];
    ellBase = [sqrt(eigVal(2,2))*cos(angles); sqrt(eigVal(1,1))*sin(angles)]';
    ellRot = eigVec(:,2:-1:1) * ellBase';   % invert order of eigenvectors as these are the opposite way round from in R
    ellRot(1,:) = ellRot(1,:) + ctr(1);
    ellRot(2,:) = ellRot(2,:) + ctr(2);
    
    totaldistfrom0 = sqrt(ellRot(1,:).^2 + ellRot(2,:).^2);
    [n,i] = min(totaldistfrom0);
    [f,j] = max(totaldistfrom0);
    
    output.lowerCI = n;
    output.upperCI = f;
    if quantiles==68
        output.lowerCI = ((output.lowerCI - output.meanamp)/1.96) + output.meanamp;
        output.upperCI = ((output.upperCI - output.meanamp)/1.96) + output.meanamp;
    end
    
end


if (method(1)=='a')
    % calculate the mean and the error bars using absolute amplitude values (incoherent averaging)
    adata = abs(compdata);
    output.meanamp = mean(adata); % mean absolute amplitude
    se = std(adata)/sqrt(length(adata));
    if (quantiles==95)
        ci = 1.96*se;
    end
    if (quantiles==68)
        ci = se;
    end
    output.lowerCI = max(output.meanamp - ci,0);  % lower bound at 0
    output.upperCI = output.meanamp + ci;
end


end