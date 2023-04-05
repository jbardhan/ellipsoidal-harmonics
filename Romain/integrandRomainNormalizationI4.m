function fcnvals = integrandRomainNormalizationI4(lambda, n, p, a, ...
						  b, c)
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);
lame = calcLame(lambda, n, p, a, b, c);
lame = lame';
fcnvals = lambda.^2 .* lame.^2 ./ (sqrt((h^2-lambda.^2).*(k^2- ...
						  lambda.^2)));
