function fcnvals = integrandForNormalization(mu,nu,n,m,a,b,c)
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

% from Romain, eq 20
weightingFcn = (mu.^2-nu.^2) ./ (sqrt((mu.^2-h^2).*(k^2-mu.^2)).*sqrt((h^2-nu.^2)*(k^2-nu.^2)));
Emu = calcLame(mu,n,m,a,b,c);
Enu = calcLame(nu,n,m,a,b,c);

fcnvals = (Emu.*Enu).^2.*weightingFcn';


