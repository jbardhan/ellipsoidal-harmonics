complex_radius = 13.4;
ligand_radius  = 2.6;
ligand_center  = complex_radius - ligand_radius;
nc_ligand = 2;  % number of ligand charges
nc_receptor = 0; % number of receptor charges
min_distance = 1.4; % minimum distance between charges;
epsIn  = 2;
epsOut = 80;

qLigandPair = [1; -1];
qLigandSelf = [1; 0];

xyzLigandPair = [1.4 0.0 0.0; -1.4 0 0;];

qComplex = qLigandPair;
xyzComplex = xyzLigandPair;

ligandPqrData = struct('q', qLigandPair, 'xyz', xyzLigandPair);
complexPqrData = struct('q', qComplex, 'xyz', xyzComplex);


Nmax = 20; % order of expansion (0 = monopole, etc)


[Le, Le_cfa, Le_p, Le_s, Le_c] = doAnalytical(ligand_radius, ...
													 epsIn, epsOut,...
													 ligandPqrData,...
													 Nmax);
Lexact=real(Le); [V,D] = eig(Lexact);
Lexact_cfa=real(Le_cfa); [Vcfa,Dcfa] = eig(Lexact_cfa);
Lexact_p=real(Le_p); [Vp,Dp] = eig(Lexact_p);

[Ce, Ce_cfa, Ce_p, Ce_s, Ce_c] = doAnalytical(complex_radius, ...
													 epsIn, epsOut,...
													 complexPqrData,...
													 Nmax);
Cexact=real(Ce); 
Cexact_cfa=real(Ce_cfa); 
Cexact_p=real(Ce_p); 

L = Cexact-Lexact;
Epair = 0.5 * 332.112 * qLigandPair' * L * qLigandPair
Eself = 0.5 * 332.112 * qLigandSelf' * L * qLigandSelf

