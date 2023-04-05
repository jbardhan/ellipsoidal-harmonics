function cartpnt = convertToCart(pnt)

cartpnt = zeros(size(pnt,1),3);
for i = 1:size(pnt,1)
  rho = pnt(i,1); theta = pnt(i,2); phi = pnt(i,3);
  x = rho * cos(theta) * sin(phi);
  y = rho * sin(theta) * sin(phi);
  z = rho * cos(phi);
  cartpnt(i,:) = [x y z];
end

