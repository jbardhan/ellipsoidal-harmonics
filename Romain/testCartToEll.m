function recap = testCartToEll(a,b,c,pnt,lambda)

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

recap_pntsquared(1) = prod(lambda.^2)/(h^2*k^2);
recap_pntsquared(2) = (lambda(1)^2-h^2)*(lambda(2)^2-h^2)*(h^2-lambda(3)^2)/h^2/(k^2-h^2);
recap_pntsquared(3) = (lambda(1)^2-k^2)*(k^2-lambda(2)^2)*(k^2-lambda(3)^2)/k^2/(k^2-h^2);

recap = sqrt(recap_pntsquared);