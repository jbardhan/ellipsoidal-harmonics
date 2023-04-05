function val = int22(tau, a, b, c)
k = sqrt((a^2-b^2)/(a^2-c^2);
kprime = sqrt(1-k^2);
t = c^2/tau^2 - c^2;
L22 = sqrt(1-t^2)*sqrt(1-k^2*t^2);
val = 1.0/L22^2/sqrt(((a^2-c^2)*tau^2+c^2)*((b^2-c^2)*tau^2+c^2));