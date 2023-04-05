function val = int01(tau, a, b, c)
x = c^2/tau^2 - c^2;
L01 = x;
val = 1.0/L01^2/sqrt(((a^2-c^2)*tau^2+c^2)*((b^2-c^2)*tau^2+c^2));