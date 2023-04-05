function lambda = newce(A, B, C, pnt)

h = sqrt(A^2-B^2)
k = sqrt(A^2-C^2)

mu2nu2 = (k*h*pnt(1)/A)^2;

yterm = pnt(2)^2*(h^2*(k^2-h^2)/(A^2-h^2));
zterm = pnt(3)^2*(k^2*(k^2-h^2)/(A^2-k^2));

r = roots([h^2 0 -(mu2nu2+yterm+h^4) 0 mu2nu2*h^2]);

candidates = r(find(r>h));
lambda(1) = A;
lambda(2) = candidates;
lambda(3) = sqrt(mu2nu2/lambda(2)^2);
