pqrFile = '../geometry/tripeptide/tripeptide.pqr';
meshFile = '../geometry/tripeptide/trip1';

epsIn  = 4;
epsOut = 80;

[geometry, dG_cfa, dG_p, dG_lb] = doBIBEEapprox(meshFile, pqrFile, ...
																epsIn, epsOut);

geometry.A = calcBEMSystemMatrix(geometry.meshData, epsIn, epsOut);
Lbem = geometry.C * (geometry.A\geometry.B);
Lbem_cfa = geometry.C * (geometry.A_cfa\geometry.B);
Lbem_p   = geometry.C * (geometry.A_p\geometry.B);
Lbem_mod = calcModBIBEE(geometry);
Lbem_interp = calcInterpBIBEE(geometry);

[Vbem,Dbem]=eig(Lbem); dbem=diag(Dbem);
[Vcfa,Dcfa]=eig(Lbem_cfa); dcfa=diag(Dcfa);
[Vp,Dp]=eig(Lbem_p); dp=diag(Dp);
[Vm,Dm]=eig(Lbem_mod); dm = diag(Dm);
[Vi,Di]=eig(Lbem_interp); di = diag(Di);