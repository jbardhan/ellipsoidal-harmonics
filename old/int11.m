function val = int11(tau, a, b, c)
t = c^2./tau.^2 - c^2;
L11 = sqrt(1-t.^2);
val = 1./(L11.^2.*(sqrt(((a^2-c^2)*tau.^2+c^2).*((b^2-c^2)*tau.^2+ ...
																 c^2))));
