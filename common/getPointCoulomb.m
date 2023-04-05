function [B,C] = getPointCoulomb(pqrData, points, normals)

numcharges = length(pqrData.q);
numpoints  = size(points,1);
C = zeros(numpoints,1);
B = C;

for i=1:numpoints
  for j=1:numcharges
	 rvec = points(i,:) - pqrData.xyz(j,:);
	 r = norm(rvec);
	 if r < 1e-10
		G = 0; dGdn = 0;
	 else
		G = 1.0/4/pi/r;
		dGdn = -dot(rvec,normals(i,:))/4/pi/r^3;
	 end
	 C(i,j) = G;
	 B(i,j) = dGdn;
  end
end
