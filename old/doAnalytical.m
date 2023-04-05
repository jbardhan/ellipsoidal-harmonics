function [L, Lcfa, Lp, Ls, Lc, Li] = doAnalytical(b, epsIn, epsOut, pqrData, Nmax)

% b is the sphere radius, in Angstroms
% pqrData contains an M-length vector q of charge values and an
% M-by-3 matrix xyz of their locations
% epsIn and epsOut are the dielectric constants
% Nmax is the maximum multipole order to use
% L is the actual Hessian
% Lbib is the BIBEE/CFA Hessian
% (to get kcal/mol energies, you need to multiply by 332.112
% outside this function)
M = length(pqrData.q);
for i=1:M
  newq = zeros(M,1);
  newq(i) = 1.0;

  Enm = computeEnm(b, epsIn, pqrData, newq, Nmax);

  [Bnm_exact,Snm,Snm2] = computeBnm_exact(b, Enm, epsIn, epsOut, Nmax);
  L(:,i) = computePot(Bnm_exact, pqrData, Nmax);

  Bnm_bibee = computeBnm_bibee(b, Enm, epsIn, epsOut, Nmax);
  Lcfa(:,i) = computePot(Bnm_bibee, pqrData, Nmax);
  
  Bnm_P = computeBnm_bibeeP(b, Enm, epsIn, epsOut, Nmax);
  Lp(:,i) = computePot(Bnm_P, pqrData, Nmax);

  Bnm_S = computeBnm_bibeeS(b, Enm, epsIn, epsOut, Nmax);
  Ls(:,i) = computePot(Bnm_S, pqrData, Nmax);

  Bnm_I = computeBnm_interp(b, Enm, epsIn, epsOut, Nmax);
  Li(:,i) = computePot(Bnm_I, pqrData, Nmax);

  Bnm_C = Bnm_P;
  Bnm_C(1,:) = Bnm_bibee(1,:);
  Lc(:,i) = computePot(Bnm_C, pqrData, Nmax);

  %  keyboard
end
