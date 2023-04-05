function [V,D,M,d,f,g] = doL(n, a,b,c)

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

alpha = h^2;
beta = k^2 - h^2;
gamma =alpha - beta;

r = floor(n/2);

iVec = 0:n-r-1; 
iVecShort = 0:n-r-2;

if ~mod(n,2)  % i.e., if n is even   Ritter type 6
  f = -(2*r-2*(iVecShort+1)).*(2*r+2*(iVecShort+1)+1)*alpha;
  d = (2*r*(2*r+1)-(2*iVec+1).^2)*alpha + beta*(2*iVec+2).^2;
  g = -beta*(2*iVecShort+2).*(2*iVecShort+3);
else  % if n is odd                  Ritter type 3
iVec = 0:n-r-1; 
iVecShort = 0:n-r-2;

  f = -alpha*(2*r-2*(iVecShort+1)+2).*(2*r+2*(iVecShort+1)+1);
  d = alpha*(2*r+1)*(2*r+2)-gamma*(2*iVec+1).^2;
  g = -beta*(2*iVecShort+2).*(2*iVecShort+3);
end
M = diag(d)+diag(f,-1)+diag(g,1);
[V,D]=eig(M);
