function E = calcLameProduct(a,b,c,srcEll,nmax)

E = zeros(nmax^2,1);
index = 1;
for n=0:nmax
  for p=0:2*n
    E_lambda(index)= calcLame(srcEll(1), n, p, a, b, c,srcEll(2),srcEll(3));
    E_mu(index)    = calcLame(srcEll(2), n, p, a, b, c,srcEll(2),srcEll(3));
    E_nu(index)    = calcLame(srcEll(3), n, p, a, b, c,srcEll(2),srcEll(3));
    E(index)  = E_lambda(index) * E_mu(index) * E_nu(index);
    index = index+1;
  end
end
