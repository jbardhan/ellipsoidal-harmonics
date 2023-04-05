pqrFile = '../geometry/sphere/sphere.pqr';
meshFile = '../geometry/sphere/sphere_r5';
sphere_radius = 5;
epsIn  = 4;
epsOut = 80;
Nmax = 30; % order of expansion (0 = monopole, etc)
nc = 50;
qnew = ones(nc,1);
xyznew=generatePoints(nc,pqrFile);
geometry.pqrData.q = qnew;
geometry.pqrData.xyz = xyznew;
[Le, Le_cfa, Le_p, Le_s, Le_c] = doAnalytical(sphere_radius, ...
													 epsIn, epsOut,...
													 geometry.pqrData,...
													 Nmax);

Lexact=real(Le); [V,D] = eig(Lexact);
Lexact_cfa=real(Le_cfa); [Vcfa,Dcfa] = eig(Lexact_cfa);
Lexact_p=real(Le_p); [Vp,Dp] = eig(Lexact_p);
Lexact_s=real(Le_s); [Vs,Ds] = eig(Lexact_s);
Lexact_c=real(Le_c); [Vc,Dc] = eig(Lexact_c);

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