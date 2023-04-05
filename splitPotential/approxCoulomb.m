function [phi, qmod, coefficients] = approxCoulomb(a,b,c,src, dest, q, order)

srcEll = approxCartToEll(a,b,c,src);
destEll = approxCartToEll(a,b,c,dest);

qmod = q;
phi  = 0;
index = 1;
gamma = computeRomainNormalizationConstants(order, a,b,c);
I = computeExteriorIntegral(destEll(1),order,a,b,c);
for n=0:order
  for p=0:2*n
    E_lambda(index)= calcLame(srcEll(1), n, p, a, b, c,srcEll(2),srcEll(3));
    E_mu(index)    = calcLame(srcEll(2), n, p, a, b, c,srcEll(2),srcEll(3));
    E_nu(index)    = calcLame(srcEll(3), n, p, a, b, c,srcEll(2),srcEll(3));
    E(index)  = E_lambda(index) * E_mu(index) * E_nu(index);

    Edest_lambda(index)= calcLame(destEll(1), n, p, a, b, c,destEll(2),destEll(3));
    Edest_mu(index)    = calcLame(destEll(2), n, p, a, b, c,destEll(2),destEll(3));
    Edest_nu(index)    = calcLame(destEll(3), n, p, a, b, c,destEll(2),destEll(3));
    Edest(index)  = Edest_lambda(index) * Edest_mu(index) * Edest_nu(index);

    coefficients(index) = ((4 * pi)/(2*n+1))*(1./gamma(index)) * E(index);
    qmod = qmod - coefficients(index)*E(index);

    F(index) = (2*n+1)*Edest(index)*I(index);
    phi  = phi + coefficients(index) * F(index);
    index = index+1;
  end
end

%keyboard