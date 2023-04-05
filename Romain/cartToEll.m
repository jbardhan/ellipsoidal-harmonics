function lambda = cartToEll(a,b,c,pnt)

lambdaApprox = approxCartToEll(a, b, c, pnt);

lambda = lambdaApprox;

return;

for i=1:3
  vals = [lambdaApprox(i)^2 
			 abs(lambdaApprox(i)^2-h^2)
			 abs(lambdaApprox(i)^2-k^2)];
  [srted,I]=sort(vals,'ascend');
  
  if I(1)==1
	 polycoeffs = [1
					  ];
  elseif I(1)==2
	 polycoeffs = [1
					  ];
  else % I(1)==3
	 polycoeffs = [1
					  ];
  end
  
end
