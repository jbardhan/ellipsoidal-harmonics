function F = computeDestExpansionMatrix(a,b,c,pts, order)

numpts = size(pts,1);
Nvector = [];
for n=0:order
  Nvector = [Nvector; ones(2*n+1,1)*(2*n+1)];
end

F = zeros(numpts, (order+1)^2);
for ptIndex = 1:numpts
  destEll = approxCartToEll(a,b,c,pts(ptIndex,:));
  I = transpose(computeExteriorIntegral(destEll(1), order, a,b,c));
  lameProducts = calcLameProduct(a,b,c,destEll,order);
  assembledProjections = Nvector .* lameProducts .* I;
  F(ptIndex,:) = transpose(assembledProjections);
end
