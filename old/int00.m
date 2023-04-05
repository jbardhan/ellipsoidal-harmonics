function val = int00(tau,a,b,c)
x = c^2./tau.^2 - c^2;
L00 = ones(1,length(tau));
val = 1./(L00.^2.*(sqrt(((a^2-c^2)*tau.^2+c^2).*((b^2-c^2)*tau.^2+ ...
																 c^2))));
