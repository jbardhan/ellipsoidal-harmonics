function [L,Lprime] = getLambdaD(a,b,c)

sumSquare = a^2+b^2+c^2;
sqrtSum = sqrt((a^4-b^2*c^2)+(b^4-a^2*c^2)+(c^4-a^2*b^2));
L      = sumSquare/3.0 + sqrtSum/3.0;
Lprime = sumSquare/3.0 - sqrtSum/3.0;