function val = int21(tau, a, b, c)
k = sqrt((a^2-b^2)/(a^2-c^2);
kprime = sqrt(1-k^2);
t = c^2/tau^2 - c^2;
L21 = t^2 - (1+k^2-sqrt(1-k^2*kprime^2))/3/k^2;
val = 1.0/L21^2/sqrt(((a^2-c^2)*tau^2+c^2)*((b^2-c^2)*tau^2+c^2));