function [L, Lprime, K,lambda] = calcLame(a,b,c, Nmax)

if Nmax > 3
  L = 0;
  K = 0;
  error('Cannot compute Lame functions for higher order than N=3!\n');
end

if ((a < b) || (a < c) || (b < c))
  L = 0;
  K = 0;
  error('Must give me an ellipsoid with a>b>c!\n');
end

k = sqrt((a^2 - b^2)/(a^2 - c^2));
kprime = sqrt(1-k^2);

L = zeros(2*Nmax+1);
Lprime = L;
K = L;

% stored in order,term L(n,m) so m ranges from 0 to 2*n
t = 0;
L(0+1,0+1) = 1;  % L_0^0 =1
L(1+1,0+1) = t;  % L_1^0 = t
L(1+1,1+1) = sqrt(1-t^2);  % L_1^1
L(1+1,2+1) = sqrt(1-k^2*t^2); % L_1^2
L(2+1,0+1) = t^2-(1+k^2+sqrt(1-k^2*kprime^2))/3.0/k^2; % L_2^0
L(2+1,1+1) = t^2-(1+k^2-sqrt(1-k^2*kprime^2))/3.0/k^2; % L_2^1
L(2+1,2+1) = sqrt(1-t^2)*sqrt(1-k^2*t^2); % L_2^2
L(2+1,3+1) = t*sqrt(1-t^2); % L_2^3
L(2+1,4+1) = t*sqrt(1-k^2*t^2); %L_2^4;

Lprime(0+1,0+1) = 0;  % Lprime_0^0 =1
Lprime(1+1,0+1) = 1;  % Lprime_1^0 = t
Lprime(1+1,1+1) = -t/sqrt(1-t^2);  % Lprime_1^1
Lprime(1+1,2+1) = -k^2*t/sqrt(1-k^2*t^2); % Lprime_1^2
Lprime(2+1,0+1) = 2*t; % Lprime_2^0
Lprime(2+1,1+1) = 2*t; % Lprime_2^1
Lprime(2+1,2+1) = -t*(sqrt(1-k^2*t^2)/sqrt(1-t^2)+k^2*sqrt(1-t^2)/sqrt(1-k^2*t^2)); % Lprime_2^2
Lprime(2+1,3+1) = (sqrt(1-t^2) + -t^2/sqrt(1-t^2)); % Lprime_2^3
Lprime(2+1,4+1) = (sqrt(1-k^2*t^2) + -k^2*t^2/sqrt(1-k^2*t^2)); %Lprime_2^4;

K(0+1,0+1) = (2*0+1) * L(0+1,0+1) * 2*c*quad(@(tau)int00(tau,a,b,c),0,1);
K(1+1,0+1) = (2*1+1) * L(1+1,0+1) * 2*c*quad(@(tau)int10(tau,a,b, ...
																  c),0,1);
K(1+1,1+1) = (2*1+1) * L(1+1,1+1) * 2*c*quad(@(tau)int11(tau,a,b,c),0,1);
K(1+1,2+1) = (2*1+1) * L(1+1,2+1) * 2*c*quad(@(tau)int12(tau,a,b, ...
																  c),0,1);
if 0
K(2+1,0+1) = (2*2+1) * L(2+1,0+1) * 2*c*quad(@(tau)int20(tau,a,b,c),0,1);
K(2+1,1+1) = (2*2+1) * L(2+1,1+1) * 2*c*quad(@(tau)int21(tau,a,b, ...
																  c),0,1);
K(2+1,2+1) = (2*2+1) * L(2+1,2+1) * 2*c*quad(@(tau)int22(tau,a,b,c),0,1);
K(2+1,3+1) = (2*2+1) * L(2+1,3+1) * 2*c*quad(@(tau)int23(tau,a,b, ...
																  c),0,1);
K(2+1,4+1) = (2*2+1) * L(2+1,4+1) * 2*c*quad(@(tau)int24(tau,a,b,c),0,1);
end

for n=0:1
  i= n+1;
  for m=0:1*i
	 j=m+1;
	 lambda(i,j) = -1 + 2*a*b*c/(2*n+1) * Lprime(i,j)*K(i,j);
  end
end
