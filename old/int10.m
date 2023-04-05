function val = int10(tau, a, b, c)
x = c^2./tau.^2 - c^2;
L10 = x;
val = 1./(L10.^2.*(sqrt(((a^2-c^2)*tau.^2+c^2).*((b^2-c^2)*tau.^2+ ...
																 c^2))));
