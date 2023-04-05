function [A,centroids,normals,panelAreas,Dstar] = calcBEMSystemMatrix(meshData,epsIn,epsOut)

numPanels = size(meshData.face,1);
panelInts = zeros(numPanels, numPanels);

for i=1:numPanels
  panelVerts = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  [panelAreas(i),centroids(i,:),normals(i,:)] = calcp(panelVerts);
end

for i=1:numPanels
  panelVerts = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  [j1,j2,j3,j4,panelInts(i,:)] = calcp(panelVerts,centroids);
end
% now panelInts(i,j) = potential induced by dipole layer on panel i 
% at the centroid of panel j

% scaling all the columns by the appropriate panel areas gives us the
% qualocation-discretized operator (see Tausch et al, 2001, Altman et
% al 2005, Bardhan 2009, Bardhan, Eisenberg, and Gillespie 2009).  the
% negative sign shows up because we are actually computing the
% electric field due to a monopole (thus taking the gradient with
% respect to the OTHER position).

A = -panelInts' * diag(panelAreas);
A = A/epsIn;
% now fix up the diagonals with the discretized identity operator
for i=1:numPanels
  A(i,i) = panelAreas(i)* 2*pi*(epsIn+epsOut)/(epsIn*(epsIn-epsOut));
end
