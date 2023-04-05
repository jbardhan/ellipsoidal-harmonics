function gamma = computeNormalizationConstants(n,a,b,c)
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

muMin = h*1.000001; muMax = 0.999999*k;
nuMin =-h*0.999999; nuMax = 0.999999*h;
index = 1;
for i=0:n
  for j=0:2*i
	 fprintf('%d: %d,%d\n',n, i,j);
	 gamma(index) = dblquad(@(mu, ...
				  nu)integrandForNormalization(mu,nu,i,j,a,b,c),muMin,muMax,nuMin,nuMax);
	 index = index+1;
  end
end

gamma = 4*gamma; % we are using the Dassios normalization constants
                 % (B16-B20)
		 