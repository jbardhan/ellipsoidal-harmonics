function [lambdaV, lambdaK] = computeLaplaceEigenvalues(radius, n)

highestOrderMultipole = ceil(sqrt(n));
lambdaV = zeros(highestOrderMultipole,1);
lambdaK = lambdaV;

for i=1:highestOrderMultipole
  lambdaV(1+(i-1)^2:i^2) = radius/(2*(i-1)+1);
  lambdaK(1+(i-1)^2:i^2) = -1./(2.*(2*(i-1)+1.0));
end


lambdaV = lambdaV(1:n);
lambdaK = lambdaK(1:n);

