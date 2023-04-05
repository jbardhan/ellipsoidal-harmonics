function [M,Mderiv] = calcM(lambda, n, a, b, c, signm, signn)

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

r = floor(n/2);

% these are not actually used for K, but included for completeness
% here.  signCKT = clever Knepley transform
if nargin < 7
  signm = 1;
  signn = 1;
end
signCKT1 = sign(lambda) * signm; % sign of mu 
signCKT2 = sign(lambda) * signn; % sign of nu

[V,D,M,d,f,g]=doM(n, a,b,c);

if size(V,1) < 1
    M = [];
    return
end

for i=1:size(V,2)
    bvec = V(:,i);

    % the following normalizes so that the highest power of lambda
    % has coefficient unity. see Romain top p242; the same
    % normalization is used by Dassios, and it is important to
    % obtain their results D9,D10 and their normalization constants
    % gamma in B16-B20.
    bvec = bvec/(bvec(end)*(-1/h^2)^(size(V,1)-1));

    P(:,i) = zeros(length(lambda),1);
    Pderiv(:,i) = zeros(length(lambda),1);
    psi(:,i) = lambda.^(1-n+2*r).*signCKT2.*sqrt(abs(lambda.^2-k^2));
    psideriv(:,i) = lambda.^(-n+2*r).*signCKT2.*sqrt(abs(lambda.^2-k^2)) ...
	.* ((1-n+2*r)+(lambda.^2./(abs(lambda.^2-k^2))));
    for j=1:size(V,1)
        P(:,i) = P(:,i)+bvec(j)*(1-lambda.^2/h^2).^(j-1);
	Pderiv(:,i)=Pderiv(:,i)+bvec(j)*(1-lambda.^2/h^2).^(j-2).*(-2*(j-1)*lambda/h^2);
    end
end
%keyboard
M = P.*psi;
Mderiv = P.*psideriv + Pderiv.*psi;