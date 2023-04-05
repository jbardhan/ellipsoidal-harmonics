function [E1,E2,E] = computeSourceExpansionMatrix(a,b,c,pts,o1,o2)

numCharges = size(pts,1);
E = zeros((o2+1)^2,numCharges);

gamma = transpose(computeRomainNormalizationConstants(o2,a,b,c));

scaleNvector = [];
for n=0:o2
  scaleNvector = [scaleNvector; 4*pi*ones(2*n+1,1)/(2*n+1)];
end

index =1;
for qIndex = 1:numCharges
  srcEll = approxCartToEll(a,b,c,pts(qIndex,:));
  lameProductsForSourcePoint = calcLameProduct(a,b,c,srcEll,o2);
  assembledProjections = scaleNvector .* (1./gamma) .* lameProductsForSourcePoint;
%  keyboard
  E(:,qIndex) = assembledProjections;
end

E1 = E(1:(o1+1)^2,:);
E2 = E((o1+1)^2+1:end,:);
