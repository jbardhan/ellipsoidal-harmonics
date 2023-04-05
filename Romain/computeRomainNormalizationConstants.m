function gamma = computeRomainNormalizationConstants(nmax, a, b, c)
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

index = 1;

hscale =  h;
kscale =  k;
hscale2 = h;
for n=0:nmax
  for p=0:2*n
% Romain Annex 2
mathcalI1 =quadgk(@(lambda)integrandRomainNormalizationI1(lambda,n, ...
						  p,a,b,c),hscale,kscale);
mathcalI2 =quadgk(@(lambda)integrandRomainNormalizationI2(lambda,n, ...
						  p,a,b,c),hscale,kscale);
mathcalI3 = quadgk(@(lambda)integrandRomainNormalizationI3(lambda,n, ...
						  p,a,b,c),0,hscale2);
mathcalI4 = quadgk(@(lambda)integrandRomainNormalizationI4(lambda,n, ...
						  p,a,b,c),0,hscale2);

% equation 53, Iij = I_i^j

I20 = (1/2 *h)*quadgk(@(Lambda)integrandRomainI2(Lambda,0,a,b,c),0,1-k^2/h^2);
I21 = (1/2 *h)*quadgk(@(Lambda)integrandRomainI2(Lambda,1,a,b,c),0,1-k^2/h^2);
I30 = (1/2 *h)*quadgk(@(Lambda)integrandRomainI3(Lambda,0,a,b,c),0,1);
I31 = (1/2 *h)*quadgk(@(Lambda)integrandRomainI3(Lambda,1,a,b,c),0,1);

matrix = [I20 I21; I30 I31];
sol1 = matrix \ [mathcalI1; mathcalI3];
alpha = sol1(1);
beta  = sol1(2);

sol2 = matrix \ [mathcalI2; mathcalI4];
A = sol2(1);
B = sol2(2);

gamma(index) = (pi/2)*(alpha*B-beta*A); % Romain eq 54
index = index+1;

  end
end


gamma = 8 * gamma; % now it is normalized according to Dassios
% see discussion in Deng paper, the different uses are related
% to integrating over the whole ellipsoid or just one octant.

