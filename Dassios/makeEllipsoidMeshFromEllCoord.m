function [tri,x,y,z,areas,centroids,normals,ellipseGrid]= ...
	 makeEllipsoidMeshFromEllCoord(xc, yc,zc,xr,yr,zr,Np)

[h,Lambda,LambdaPrime,Gamma,alpha]=doDassiosEllipsoid([xr yr zr]);
rho = alpha(1);
[xs,ys,zs,MU,NU]=partialEllipsoidFromEllipsoidCoord(alpha,h,Np);
ellipseGrid = struct('rho',rho,'MU',MU,'NU',NU,'xs',xs,'ys',ys,'zs',zs,'alpha',alpha,'Gamma',Gamma,'h',h,'Lambda',Lambda,'LambdaPrime',LambdaPrime);
np = prod(size(xs));
[t1,x1,y1,z1]=maketri(xs,ys,zs);
[t2,x2,y2,z2]=makefliptri(xs,-ys,zs);
[t3,x3,y3,z3]=maketri(xs,-ys,-zs);
[t4,x4,y4,z4]=makefliptri(xs,ys,-zs);
%t=t1; x=x1;y=y1;z=z1;
t = [t1;t2+np;t3+2*np;t4+3*np];
x = [x1;x2;x3;x4];
y = [y1;y2;y3;y4];
z = [z1;z2;z3;z4];

[tri,areas,centroids,normals]=prunemesh(t,x,y,z);
