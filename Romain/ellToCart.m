function pnt = ellToCart(a,b,c,lambda,mu,nu)

if nargin < 6
    mu = lambda(2);
    nu = lambda(3);
    lambda = lambda(1);
end
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

signmu = sign(mu);
signnu = sign(nu);

n = 1; r = floor(n/2);
Kpsi_n = calcK(lambda, n, a, b, c,signmu,signnu)*calcK(mu, n, a, b, c,signmu,signnu)*calcK(nu, n, a, b, c,signmu,signnu);
Lpsi_n = calcL(lambda, n, a, b, c,signmu,signnu)*calcL(mu, n, a, b, c,signmu,signnu)*calcL(nu, n, a, b, c,signmu,signnu);
Mpsi_n = calcM(lambda, n, a, b, c,signmu,signnu)*calcM(mu, n, a, b, c,signmu,signnu)*calcM(nu, n, a, b, c,signmu,signnu);

pnt(1) = (lambda^2*mu^2*nu^2)/h^2/k^2;
pnt(2) = (lambda^2-h^2)*(mu^2-h^2)*(h^2-nu^2)/h^2/(k^2-h^2);
pnt(3) = (lambda^2-k^2)*(k^2-mu^2)*(k^2-nu^2)/k^2/(k^2-h^2);
pnt = sqrt(pnt);

if Kpsi_n < 0
  pnt(1) = -pnt(1);
end
if Lpsi_n < 0
  pnt(2) = -pnt(2);
end
if Mpsi_n < 0
  pnt(3) = -pnt(3);
end
