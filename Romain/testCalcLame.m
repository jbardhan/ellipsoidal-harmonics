testOrder=4;
a=3;
b=2;
c=1;
doLambda = 1;
numLambda = 50;
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

n = testOrder;
r = floor(n/2);

if doLambda == 1
    lambda = sqrt(k^2:(9*k^2)/numLambda:10*k^2);
elseif doLambda == 2
    lambda = sqrt(h^2:(k^2-h^2)/numLambda:k^2);
elseif doLambda == 3
    lambda = sqrt(0:h^2/numLambda:h^2);
end
lambda = lambda';

K = calcK(lambda, testOrder, a, b, c);
L = calcL(lambda, testOrder, a, b, c);
M = calcM(lambda, testOrder, a, b, c);
N = calcN(lambda, testOrder, a, b, c);
