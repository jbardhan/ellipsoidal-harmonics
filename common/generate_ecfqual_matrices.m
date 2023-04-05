function [A,B,C] = generate_ecfqual_matrices(surfaceToSurfaceOperators, ...
					     chargeToSurfaceOperators, ...
					     surfaceToChargeOperators, ...
					     areas, E_0,E_inf,E_in, E_out)
% note that the matrices A and B are missing a factor of (1/E_in)
% compared to Altman05 eqs 18-20... it will after all cancel out.

A = surfaceToSurfaceOperators.K' * diag(areas);
A = A + 0.5 * ((E_in+E_out)/(E_in-E_out)) *diag(areas);
B = -surfaceToChargeOperators.dlpToCharges'; 
C = surfaceToChargeOperators.slpToCharges/E_in;
