% this tests BEM calculation of the exact shape. it does NOT give the
% exact ellipsoidal result, e.g. from Xue and Deng (Phys. Rev. E,
% v.83:056709 (2011)).

addpath('../common','../elliptic_package','../Dassios');
test = 1;

pqr = readpqr('../geometry/ellipsoid_points/monopole.pqr');

loadconstants;
R = 1;
E_protein = 4;
E_water = 80;
kcalmol = 332.112;

if test == 1
  alpha = [7.354276 5.2282265 3.489720];
else
  alpha = [5.815856 3.103770 2.325190];
end

Ntest = [10];

for i=1:length(Ntest)
[tri,x,y,z,areas,centroids,normals]=makeEllipsoidMesh(0,0,0, ...
						  alpha(1),alpha(2), ...
						  alpha(3), ...
						  Ntest(i));
meshX = x(tri)'; meshY = y(tri)'; meshZ = z(tri)';
meshData = struct('face',tri,'X',meshX,'Y',meshY,'Z',meshZ);
Npanels(i) = length(areas);

[A] = collocation_mesh(meshData,centroids,normals,areas);
[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas, pqr);
[Aqual,Bqual,Cqual]=generate_ecfqual_matrices(A,B,C,areas,E_0, ...
					      E_inf,E_protein,E_water);

[A_cfa, A_P, A_M] = generateEllipsoidBIBEE(A,B,C,areas,E_protein,E_water,Aqual);

E_bem(i) = 0.5 * kcalmol * 4 * pi * (pqr.q' * Cqual * (Aqual \ (Bqual* pqr.q)))
E_P(i) = 0.5 * kcalmol * 4 * pi * (pqr.q' * Cqual * (A_P \ (Bqual * ...
						  pqr.q)));
E_cfa(i)= 0.5 * kcalmol * 4 * pi * (pqr.q' * Cqual * (A_cfa \ (Bqual * ...
						  pqr.q)));
E_M(i) = 0.5 * kcalmol * 4 * pi * (pqr.q' * Cqual * (A_M \ (Bqual ...
						  * pqr.q)));

end


