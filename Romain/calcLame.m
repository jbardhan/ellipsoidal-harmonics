function [lame,lameDeriv,lameMatrix,lameMatrixDeriv] = calcLame(lambda, ...
						  n, p, a, b, c, ...
						  mu, nu)

if nargin < 8
  signmu = 1;
  signnu = 1;
else
  signmu = sign(mu);
  signnu = sign(nu);
end
r = floor(n/2);
if size(lambda,2) > size(lambda,1)
  lambda = lambda';
end

if ((p >= 0) && (p < r + 1))
%  fprintf('order %d: returning %d entry of %d K functions\n',n,p+1, r+1);
  [lameMatrix,lameMatrixDeriv] = calcK(lambda, n, a, b, c,signmu,signnu);
  lame = lameMatrix(:,p+1);
  lameDeriv = lameMatrixDeriv(:,p+1);
elseif ((p >= r+1) && (p < n + 1))
  [lameMatrix,lameMatrixDeriv] = calcL(lambda, n, a, b, c,signmu,signnu);
  lame = lameMatrix(:,p-(r+1)+1);
  lameDeriv = lameMatrixDeriv(:,p-(r+1)+1);
%  fprintf('returning %d entry of %d L functions\n',(p-(r+1)+1),n-r);
elseif ((p >= n+1) && (p < 2*n - r+ 1))
  [lameMatrix,lameMatrixDeriv] = calcM(lambda, n, a, b, c,signmu,signnu);
  lame = lameMatrix(:,p-(n+1)+1);
  lameDeriv = lameMatrixDeriv(:,p-(n+1)+1);
%  fprintf('returning %d entry of %d M functions\n',(p-(n+1)+1),n- ...
%			 r);
elseif ((p >=2*n-r +1) && (p < 2*n+1))
  [lameMatrix,lameMatrixDeriv] = calcN(lambda, n, a, b, c,signmu,signnu);
  lame = lameMatrix(:,p-(2*n-r+1)+1);
  lameDeriv = lameMatrixDeriv(:,p-(2*n-r+1)+1);
%  fprintf('returning %d entry of %d N functions\n',(p-(2*n-r+1)+1), ...
%			 r);
else
  fprintf('There is no (%d, %d) pair\n',n,p);
end
