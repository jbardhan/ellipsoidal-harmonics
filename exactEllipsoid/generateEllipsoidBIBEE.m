function [A_cfa, A_P, A_M] = generateEllipsoidBIBEE(A,B,C,areas, ...
						  epsIn,epsOut,Aq)

dielecFunc = 2.0*(epsIn-epsOut)/(epsIn+epsOut);
scaledI = 1./dielecFunc - 1.0/2.0;
scaledIplusLimit = scaledI + 1.0/2.0;
for i=1:length(areas)
p_vec(i)   = areas(i) * scaledIplusLimit;
cfa_vec(i) = areas(i) * scaledI;
end
A_P = diag(p_vec);
A_cfa = diag(cfa_vec);

I = eye(length(areas));

order = 16; %size(A.K,1);
[V,D]=eig(A.K');
Lred = diag(D); 
Lred = Lred(1:order);
Vred = V(:,1:order);
invV = inv(V);
invVred = invV(1:order,:);

A_M = (scaledIplusLimit*I + real(Vred * diag(Lred) * invVred) ) * diag(areas);
keyboard