function output = getel(input)
% getel: helper function that calculates an ellipse bounding some data points
% the outline of the ellipse is returned, sampled at 200 points
% the input should be an Nx2 matrix of observations or an N-length vector of complex values
% this function is part of the FourierStats package: https://github.com/bakerdh/FourierStats

compdata = input;
if (~isreal(input))
    temp(:,1) = real(compdata); 
    temp(:,2) = imag(compdata);
    compdata = temp;
end

A = cov(compdata);
ctr = mean(compdata);
RR = chol(A);
angles = linspace(0,2*pi,200);  % angles for ellipse
ell = 1 * [cos(angles); sin(angles)]' * RR;
ellCtr = ell + ctr;
[eigVec,eigVal] = eig(A);
eigScl = eigVec * sqrt(eigVal);
xMat = [ctr(1) + eigScl(1,:); ctr(1) - eigScl(1,:)];
yMat = [ctr(2) + eigScl(2,:); ctr(2) - eigScl(2,:)];
ellBase = [sqrt(eigVal(2,2))*cos(angles); sqrt(eigVal(1,1))*sin(angles)]';
ellRot = eigVec(:,2:-1:1) * ellBase';   % invert order of eigenvectors as these are the opposite way round from in R
output(1,:) = ellRot(1,:) + ctr(1);
output(2,:) = ellRot(2,:) + ctr(2);

end