function [A_cfa, A_p, A_lb, A_interp] = calcBIBEEmatrices(pqr, mesh, epsIn, epsOut)

numCharges = length(pqr.q);
numPanels  = size(mesh.face,1);

cfa_vec = zeros(numPanels,1);
p_vec   = zeros(numPanels,1);
lb_vec  = zeros(numPanels,1);

for i=1:numPanels
  panelVerts = [mesh.X(:,i) mesh.Y(:,i) mesh.Z(:,i)];
  area(i)    = calcp(panelVerts);
  % the below is the same exact expression from calcBEMSystemMatrix
  dielecFunc = 2.0*(epsIn-epsOut)/(epsIn+epsOut);
  p_vec(i)   = area(i) * (4*pi/epsIn) * 1./dielecFunc;
  interp_vec(i) = area(i) * (4*pi/epsIn) * (1./dielecFunc - 0.12);
  cfa_vec(i) = area(i) * (4*pi/epsIn) * (1./dielecFunc - 1/2);
  lb_vec(i)  = area(i) * (4*pi/epsIn) * (1./dielecFunc + 1/2);
  % see Bardhan, Knepley, Anitescu for the above expressions
end

A_cfa = diag(cfa_vec);
A_p   = diag(p_vec);
A_lb  = diag(lb_vec);
A_interp = diag(interp_vec);