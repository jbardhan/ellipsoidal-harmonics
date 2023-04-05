function [x,y,z,MU,NU]=partialEllipsoidFromEllipsoidCoord(alpha,h,Npoints)

mu= (h(2)+h(3))/2 + ((h(3)-h(2))/2) * cos((0:Npoints-1)*pi/(Npoints-1));
nu= (-h(3)+h(3))/2 + ((h(3)+h(3))/2) * cos((0:Npoints-1)*pi/(Npoints-1));
%linspace(h(2),h(3),Npoints);
%linspace(-h(3),h(3),Npoints);
[MU, NU] = meshgrid(mu,nu);
rho = alpha(1);
x = rho*(MU.*NU)/h(2)/h(3);
y = sqrt(rho^2-h(3)^2)*sqrt(MU.^2-h(3)^2).*sqrt(h(3)^2-NU.^2)/h(1)/h(3);
z = sqrt(rho^2-h(2)^2)*sqrt(h(2)^2-MU.^2).*sqrt(h(2)^2-NU.^2)/h(1)/h(2);

