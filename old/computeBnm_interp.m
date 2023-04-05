function [Bnm,Snm,Snm2] = computeBnm_interp(b, Enm, epsIn, epsOut, Nmax)

effectiveEigenvalue = -0.12;
epsHat = (epsIn-epsOut)/(0.5*(epsIn+epsOut));
Bnm = 0 * Enm;
for n=0:Nmax-1
  iIndex=n+1;

  Vlambda = b/(1+2*n);
  Klambda = -1/(2*(1+2*n));
  
  for m=-n:n
	 jIndex=m+n+1;

	 if n==0 % note that I hardcoded Klambda for n==0
		Snm(iIndex,jIndex) = epsHat/(1 + epsHat*-0.5) * (n+1)/b^(n+2)* ...
			 Enm(iIndex,jIndex);
	 else
		Snm(iIndex,jIndex) = epsHat/(1 + epsHat*effectiveEigenvalue) * (n+1)/b^(n+2)* ...
			 Enm(iIndex,jIndex);
	 end
  end
  
  Bnm(iIndex,1:jIndex) = Vlambda * b^(-n) * Snm(iIndex,1:jIndex) / epsIn;
end
%keyboard
