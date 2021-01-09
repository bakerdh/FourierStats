% test the pairwisemahal function using simulated data with known effect sizes


N = 10000;

dataA = 1*randn(N,2);
dataB = 1*randn(N,2);
dataB(:,1) = dataB(:,1) + 1;
grouplabels = ceil((1:(N*2))./N);
data = [dataA; dataB];
output = pairwisemahal(data,grouplabels);
output.D
