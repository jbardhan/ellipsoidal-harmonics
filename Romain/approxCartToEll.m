function lambdaReturn = approxCartToEll(A, B, C, x)
pointErrorThreshold = 1e-6;

h = sqrt(A^2-B^2);
k = sqrt(A^2-C^2);

a(1) = -(sum(x.^2) + h^2 + k^2);
a(2) = x(1)^2*(h^2+k^2) + x(2)^2*k^2 + x(3)^2*h^2 + h^2*k^2;
a(3) = -x(1)^2*h^2*k^2;

Q = (a(1)^2 - 3*a(2))/9.0;
R = (9*a(1)*a(2)-27*a(3)-2*a(1)^3)/54.0;
theta = acos(R / sqrt(Q^3));

lambdaSquared(1) = cos(theta/3);
lambdaSquared(2) = cos(theta/3 + 4*pi/3);
lambdaSquared(3) = cos(theta/3 + 2*pi/3);
lambdaSquared = 2*sqrt(Q)*lambdaSquared - a(1)/3;

n = 1;
lambdaApprox = sqrt(lambdaSquared);
if lambdaApprox(3) < 1e-10
  lambdaApprox(3) = 1e-8;
end

index = 1;
for i=0:1
  lambda = lambdaApprox(1)*(-1)^i;
  for j=0:1
    mu = lambdaApprox(2)*(-1)^j;
    signmu = sign(mu);
    for k=0:1
      nu = lambdaApprox(3)*(-1)^k;
      signnu = sign(nu);
      Kpsi_n = calcK(lambda, n, A, B, C,signmu,signnu)*calcK(mu, n, A, B, C,signmu,signnu)*calcK(nu, n, A, B, C,signmu,signnu);
      Lpsi_n = calcL(lambda, n, A, B, C,signmu,signnu)*calcL(mu, n, A, B, C,signmu,signnu)*calcL(nu, n, A, B, C,signmu,signnu);
      Mpsi_n = calcM(lambda, n, A, B, C,signmu,signnu)*calcM(mu, n, A, B, C,signmu,signnu)*calcM(nu, n, A, B, C,signmu,signnu);
      
      distances(index) = norm(x - ellToCart(A,B,C,real([lambda mu ...
		    nu])));
      index = index+1;
    end
  end
end


[foo, bestIndex] = min(distances);
i = (bitand(bestIndex-1,4)>0);
j = (bitand(bestIndex-1,2)>0);
k = (bitand(bestIndex-1,1)>0);

lambdaReturn = real([lambdaApprox(1)*(-1)^i lambdaApprox(2)*(-1)^j lambdaApprox(3)*(-1)^k]);

if pointErrorThreshold < min(distances);
  fprintf('Warning in approxCartToEll: could not find set of signs!\n')
  fprintf(' Best = %f\n', min(distances));
  keyboard
end

if pointErrorThreshold < norm(x-ellToCart(A,B,C,lambdaReturn)) 
  fprintf('Failed to give accurate transform\n');
  keyboard
end
