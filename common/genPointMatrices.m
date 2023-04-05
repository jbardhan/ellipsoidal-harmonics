function [V, K] = genPointMatrices(points, normals, weights)

N = length(weights);
V = zeros(N);
K = V;
for i=1:N
  for j=1:N
	 vecr = points(i,:)-points(j,:);
	 r = norm(vecr);
	 if i ~= j
		V(i,j) = weights(j)* 1/4/pi/r;
		K(i,j) = weights(j)* vecr*normals(j,:)'/4/pi/r/r/r;
	 end
  end
end

