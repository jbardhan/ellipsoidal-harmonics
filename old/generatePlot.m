pqrFile = '../geometry/sphere/sphere.pqr';
meshFile = '../geometry/sphere/sphere_r5';
sphere_radius = 10;
epsIn  = 1;
epsOut = 10;
Nmax = 20; % order of expansion (0 = monopole, etc)
nc = 25;
qnew = ones(nc,1);
numRuns = 1;
%xyznew=generatePoints(nc,pqrFile);
for i=1:numRuns
qrand = rand(nc,1)-0.5;
geometry.pqrData.q = qnew;
geometry.pqrData.xyz = xyznew;
[Le, Le_cfa, Le_p, Le_s, Le_c,S1,S2,Li] = doAnalytical(sphere_radius, ...
													 epsIn, epsOut,...
													 geometry.pqrData,...
													 Nmax);

Lexact=real(Le); [V,D] = eig(Lexact);
Lexact_cfa=real(Le_cfa); [Vcfa,Dcfa] = eig(Lexact_cfa);
Lexact_p=real(Le_p); [Vp,Dp] = eig(Lexact_p);
Lexact_s=real(Le_s); [Vs,Ds] = eig(Lexact_s);
Lexact_c=real(Le_c); [Vc,Dc] = eig(Lexact_c);
Lexact_i=real(Li); [Vi,Di] = eig(Lexact_i);
return;
[geometry, dG_cfa, dG_p, dG_lb] = doBIBEEapprox(meshFile, pqrFile, ...
																epsIn, epsOut);
geometry.pqrData.q = qnew;
geometry.pqrData.xyz = xyznew;

[geometry.B,geometry.C] = calcBEMmatrices(geometry.pqrData, ...
														geometry.meshData,...
														epsIn);

geometry.A = calcBEMSystemMatrix(geometry.meshData, epsIn, epsOut);
Lbem = geometry.C * (geometry.A\geometry.B);
Lbem_cfa = geometry.C * (geometry.A_cfa\geometry.B);
Lbem_p   = geometry.C * (geometry.A_p\geometry.B);
Lbem_mod = calcModBIBEE(geometry);
kcalConv = 332.112;
Eexact(i) = kcalConv*0.5 * qrand' * Lexact * qrand;
Eexact_cfa(i) = kcalConv*0.5 * qrand' * Lexact_cfa * qrand;
Eexact_p(i)   = kcalConv*0.5 * qrand' * Lexact_p * qrand;
Eexact_c(i)   = kcalConv*0.5 * qrand' * Lexact_c * qrand;

Ebem(i)   = kcalConv*0.5 * qrand' * Lbem * qrand;
Ebem_cfa(i) = kcalConv*0.5 * qrand' * Lbem_cfa * qrand;
Ebem_p(i)   = kcalConv*0.5 * qrand' * Lbem_p * qrand;
Ebem_c(i)   = kcalConv*0.5 * qrand' * Lbem_mod * qrand;
i
end
