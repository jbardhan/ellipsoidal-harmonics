addpath('../Romain');
addpath('../elliptic_package');
addpath('../Dassios');

a = 12; 
b = 10;
c = 9;
bigRadius = 20;

chargeSpacing = 2;
[gridX,gridY,gridZ] = meshgrid(-a:chargeSpacing:a+chargeSpacing,-b:chargeSpacing:b+chargeSpacing,-c:chargeSpacing:c+chargeSpacing);
SIZ = [prod(size(gridX)) 1];
origgridpts = [reshape(gridX,SIZ) reshape(gridY,SIZ) reshape(gridZ,SIZ)];
accepted = [];
for i=1:prod(size(gridX))
  if (1-0.1> sum(origgridpts(i,:).^2 ./ [a b c].^2))
  accepted = [accepted; i];
  end
end
gridpts = origgridpts(accepted,:);
numChargePoints = length(accepted);
q = 1.8 * (rand(numChargePoints,1)-0.5);

maxOrder = 3;

numFieldPoints = 50;
for i=1:numFieldPoints
  fieldPoints(i,:) = [0 0 0];
  while (1+0.1 > sum(fieldPoints(i,:).^2 ./ [a b c].^2))
    fieldPoints(i,:) = bigRadius * (rand(1,3) - 0.5);
  end
end

for i=1:numFieldPoints
  for j=1:numChargePoints
    A(i,j) = 1/norm(fieldPoints(i,:)-gridpts(j,:));
  end
end

order1 = 2; order2 = 8;
[E1,E2,E] = computeSourceExpansionMatrix(a,b,c,gridpts,order1, ...
				       order2);

F = computeDestExpansionMatrix(a,b,c,fieldPoints,order2);

return;

phiExact = 1.0/norm(fieldPoint - chargeLoc)

for i=0:maxOrder
  [phiHarmonicComponent(i+1), qMod(i+1),coeffs] = approxCoulomb(a,b,c, ...
			 			  chargeLoc, fieldPoint, q, i);
  phiApprox(i+1) = phiHarmonicComponent(i+1) + qMod(i+1)*1.0/norm(fieldPoint-chargeLoc);
end
%semilogy(abs(phiHarmonicComponent - phiExact),'--x');
