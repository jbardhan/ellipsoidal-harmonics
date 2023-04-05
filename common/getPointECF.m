function [A, B, C, K] = getPointECF(pqrData, points, normals, areas, ...
					  epsIn, epsOut)

[Balmost,Calmost] = getPointCoulomb(pqrData, points, normals);
B = -diag(areas) * 4* pi*Balmost/epsIn;
C = Calmost'/epsIn * 4*pi*diag(areas);

[V, K] = genPointMatrices(points, normals, areas);
A = K' *4*pi* diag(areas)/epsIn;
A = A + 2*pi*((epsIn+epsOut)/(epsIn*(epsIn-epsOut)))*diag(areas);

