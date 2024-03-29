function [geometryStruct, dG_cfa, dG_p, dG_lb] = doBIBEEapprox(meshfile, ...
																  pqrfile,epsIn,epsOut)

pqrData = readpqr(pqrfile);
meshData= readmesh(meshfile,1);
[B, C] = calcBEMmatrices(pqrData, meshData,epsIn); % NxN matrix not
                                             % computed automatically
[A_cfa, A_p, A_lb, A_interp] = calcBIBEEmatrices(pqrData, meshData,epsIn,epsOut);

geometryStruct = struct('pqrData', pqrData,...
								'meshData', meshData',...
								'B', B, 'C', C, ...
								'A_cfa', A_cfa, ...
								'A_p', A_p,...
								'A_lb', A_lb,...
								'A_interp', A_interp,...
								'A', 0);

[dG_cfa, dG_p, dG_lb] = computeBIBEE(geometryStruct, pqrData.q);
