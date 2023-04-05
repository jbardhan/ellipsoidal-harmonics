function val = int24(tau, a, b, c)
k = sqrt((a^2-b^2)/(a^2-c^2);
kprime = sqrt(1-k^2);
t = c^2/tau^2 - c^2;
L24 = t*sqrt(1-k^2*t^2);
val = 1.0/L24^2/sqrt(((a^2-c^2)*tau^2+c^2)*((b^2-c^2)*tau^2+c^2));