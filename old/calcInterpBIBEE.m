function L = calcInterpBIBEE(geometry)

nc = length(geometry.pqrData.q);
L = zeros(nc,nc);
np = size(geometry.A_cfa,1);
v = ones(np,1); v = v/norm(v);

for i=1:nc
  q = zeros(nc,1);
  q(i) =1;
  Bq = geometry.B*q;
  Bqmean = v * (v'*Bq);
  Bqmod  = Bq-Bqmean;

  sigmaCFA = geometry.A_cfa\Bqmean;
  sigmaInterp   = geometry.A_interp\Bqmod;
  L(:,i) = geometry.C * (sigmaCFA+sigmaInterp);
end
