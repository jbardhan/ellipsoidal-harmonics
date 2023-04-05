function Bnm = computeBnm_bibee(b, Enm, epsIn, epsOut, Nmax)

for n=0:Nmax-1
  iIndex=n+1;
  for m=-n:n
	 jIndex=m+n+1;
	 Bnm(iIndex,jIndex) = ...
		  ((epsIn-epsOut)/(epsIn*epsOut)) * (n+1)/(2*n+1)...
		  * (1/b^(2*n+1)) * Enm(iIndex,jIndex);
  end
end
