complex_radius = 36;
ligand_radius  = 16;
ligand_center  = complex_radius - ligand_radius;
nc_ligand = 10;  % number of ligand charges
nc_receptor = 10; % number of receptor charges
min_distance = 1.4; % minimum distance between charges;
epsIn  = 2;
epsOut = 80;

[qLigand, xyzLigand] = makeLigandChargeDistribution(ligand_radius, ...
																  nc_ligand, min_distance);
[qComplex, xyzComplex] = makeComplexChargeDistribution(complex_radius,...
																  nc_receptor, ...
																  min_distance, ...
																  ligand_center, ...
																  ligand_radius, ...
																  qLigand, xyzLigand);

ligandPqrData = struct('q', qLigand, 'xyz', xyzLigand);
complexPqrData = struct('q', qComplex, 'xyz', xyzComplex);


Nmax = 100; % order of expansion (0 = monopole, etc)


[Le, Le_cfa, Le_p, Le_s, Le_m, Le_i] = doAnalytical(ligand_radius, ...
													 epsIn, epsOut,...
													 ligandPqrData,...
													 Nmax);
Lexact=real(Le); [V,D] = eig(Lexact);
Lexact_cfa=real(Le_cfa); [Vcfa,Dcfa] = eig(Lexact_cfa);
Lexact_p=real(Le_p); [Vp,Dp] = eig(Lexact_p);
Lexact_m=real(Le_m); [Vm,Dm] = eig(Lexact_m);
Lexact_i=real(Le_i); [Vi,Di] = eig(Lexact_i);

[Ce, Ce_cfa, Ce_p, Ce_s, Ce_m,Ce_i] = doAnalytical(complex_radius, ...
																  epsIn, epsOut,...
																  complexPqrData,...
																  Nmax);
Cexact=real(Ce); 
Cexact_cfa=real(Ce_cfa); 
Cexact_p=real(Ce_p); 
Cexact_i=real(Ce_i); 
Cexact_m=real(Ce_m); 

[C11,C12,C21,C22] = splitMatrix(Cexact, nc_ligand, nc_receptor);
[C_cfa_11,C_cfa_12,C_cfa_21,C_cfa_22] = splitMatrix(Cexact_cfa, nc_ligand, nc_receptor);
[C_p_11,C_p_12,C_p_21,C_p_22] = splitMatrix(Cexact_p, nc_ligand, ...
														  nc_receptor);
[C_m_11,C_m_12,C_m_21,C_m_22] = splitMatrix(Cexact_m, nc_ligand, ...
														  nc_receptor);
[C_i_11,C_i_12,C_i_21,C_i_22] = splitMatrix(Cexact_i, nc_ligand, ...
														  nc_receptor);

%C11 = reaction pot matrix for ligand charges in complex state
%C12,C21 = solvation component of ligand-receptor interaction
%C22 = reaction pot matrix for receptor charges in complex state
