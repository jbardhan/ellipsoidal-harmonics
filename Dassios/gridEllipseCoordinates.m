Npoints = 50;
mu=linspace(h(2),h(3),Npoints);
nu=linspace(-h(3),h(3),Npoints);
[MU, NU] = meshgrid(mu,nu);
rho = alpha(1);
xe = rho*(MU.*NU)/h(2)/h(3);
ye = sqrt(rho^2-h(3)^2)*sqrt(MU.^2-h(3)^2).*sqrt(h(3)^2-NU.^2)/h(1)/h(3);
ze = sqrt(rho^2-h(2)^2)*sqrt(h(2)^2-MU.^2).*sqrt(h(2)^2-NU.^2)/h(1)/h(2);

