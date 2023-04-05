function [G,Dstar,panelAreas,Gp,Dp] = getG(meshData)

numpanels = size(meshData.face,1);
G = zeros(numpanels,numpanels);
D = G;

for i=1:numpanels
  panelVerts = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  [panelAreas(i),centroids(i,:),normals(i,:)] = calcp(panelVerts);
end

for i=1:numpanels
  panelVerts = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  [j1,j2,j3,G(i,:),D(i,:)] = calcp(panelVerts,centroids);
end

% galerkin like
G = G * diag(panelAreas);
Dstar = -D' * diag(panelAreas);
Dstar = Dstar + diag(panelAreas)*2*pi;
Gp = 0 * G;
Dp = 0 * Dstar;

% now do pointlike
for i=1:numpanels
  for j=1:numpanels
	 if i~=j
		rVec = centroids(i,:)-centroids(j,:);
		rLen = norm(rVec);
		Gp(i,j) = 1/rLen;
		Dp(i,j) = normals(i,:) * rVec'/ rLen^3;
	 else
		Gp(i,j) = 0;
		Dp(i,j) = 0;
	 end
  end
end

