function fcnvals = integrandRomainI20(Lambda,a,b,c);
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

i = 2, j = 0;
fcnvals = (Lambda.^j) ./ (sqrt(1-Lambda).*sqrt(h^2*Lambda+(k^2-h^2)).*sqrt(-Lambda));
