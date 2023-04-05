function [B, C] = calcBEMmatrices(pqrData,meshData,epsIn)
% this script and calcBEMSystemMatrix use calcp.m, which was
% originally written by Prof. Jacob K. White at MIT.  We are
% grateful for his permission to use it.

numPanels = size(meshData.face,1)
numCharges = length(pqrData.q);

% B maps charge values at the pqr.xyz locations to the (integral
% of the) normal electric field indcued at the boundary elements
B = zeros(numPanels, numCharges);
C = B;

for i=1:numPanels
  panelVerts = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  % each row of panelVerts is a matrix
  [j1,j2,j3,C(i,:),B(i,:)] = ...
		calcp(panelVerts, pqrData.xyz);
end
B = B/epsIn;
C = C'/epsIn;

