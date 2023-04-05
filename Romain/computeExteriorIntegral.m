function I = computeExteriorIntegral(lambda, n, a, b, c)

index = 1;
for i=0:n
  for j=0:2*i
    if  lambda > 0 
%	 fprintf('%d: %d,%d\n',n, i,j);

	 I(index) = quad(@(t)integrandForExterior(t,i,j,a,b,c),0, ...
			 1./lambda);
    else 
	 I(index) = quad(@(t)integrandForExterior(t,i,j,a,b,c),1./lambda, ...
			 0);
%	 fprintf('called negative lambda function!\n');
%	 keyboard
    end      
	 index = index+1;
  end
end
