function [V,gamma] = getSurfaceHarmonic(a,b,c,nmax,areas,mu,nu)

index = 1;
gamma = computeRomainNormalizationConstants(nmax,a,b,c);
for n=0:nmax
  for p=0:2*n
    Emu = calcLame(mu,n,p,a,b,c);
    Enu = calcLame(nu,n,p,a,b,c);
    V(:,index) = Emu.*Enu;%/sqrt(gamma(index));
    index = index+1;
    keyboard
  end
end
