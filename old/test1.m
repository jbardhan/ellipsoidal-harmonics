function [h,Lambda,LambdaPrime,Gamma,alpha] = test1(alpha)

alpha = sort(alpha,'descend');

h(1) = sqrt(alpha(2)^2-alpha(3)^2);
h(2) = sqrt(alpha(1)^2-alpha(3)^2);
h(3) = sqrt(alpha(1)^2-alpha(2)^2);
% nearly spherical gives h -> 0
% if alpha2=alpha1, we get h3 = 0: 
% if alpha2=alpha3, we get h1 = 0

sumB14 = 0;
for i=1:3
  sumB14=sumB14+(alpha(i)^4-(prod(alpha.^2)/alpha(i)^2));
end
Lambda = (1/3) * (sum(alpha.^2)+sqrt(sumB14));
LambdaPrime = (1/3) * (sum(alpha.^2)-sqrt(sumB14));


Gamma = zeros(3,6);

n = 0; % zero degree
Gamma(1+n,1) = 4*pi;

n = 1; % first degree
for i=1:3
  Gamma(1+n,i) = (4*pi/3)*prod(h.^2)/h(i)^2;
end

n = 2; % second degree
Gamma(1+n,1)= -(8*pi/5) *(Lambda-LambdaPrime)*prod(Lambda- ...
																  alpha.^2);
Gamma(1+n,2)= +(8*pi/5) *(Lambda-LambdaPrime)*prod(LambdaPrime-alpha.^2);

for i=1:3
  Gamma(1+n,2+i) = (4*pi/15)*h(i)^2*prod(h.^2);%
end


