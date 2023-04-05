% this tests BEM calculation of the exact shape. it does NOT give the
% exact ellipsoidal result, e.g. from Xue and Deng (Phys. Rev. E,
% v.83:056709 (2011)).

addpath('../common','../elliptic_package','../Dassios');

pqr = readpqr('../geometry/ion/ion.pqr');

loadconstants;
R = 1;
E_protein = 4;
E_water = 80;
kcalmol = 332.112;

alpha = [R R 2 * R];
Ntest = [3 5 10 15];

for i=1:length(Ntest)
[tri,x,y,z,areas,centroids,normals]=makeEllipsoidMesh(0,0,0, ...
																  alpha(1),alpha(2), ...
																  alpha(3),Ntest(i));
meshX = x(tri)'; meshY = y(tri)'; meshZ = z(tri)';
meshData = struct('face',tri,'X',meshX,'Y',meshY,'Z',meshZ);
Npanels(i) = length(areas);

[A] = collocation_mesh(meshData,centroids,normals,areas);
[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas, ...
													pqr);
[Aqual,Bqual,Cqual]=generate_ecfqual_matrices(A,B,C,areas,E_0, ...
															 E_inf,E_protein,E_water);

E(i) = 0.5 * kcalmol * 4 * pi * (pqr.q' * Cqual * (Aqual \ (Bqual* ...
																  pqr.q)))
end


figure;
loglog(Npanels,abs(E-analytical),'-o');
set(gca,'fontsize',16);
xlabel('Number of triangles');
ylabel('Deviation from analytical \Delta G_{solv} (kcal/mol)');

