function [tri,x,y,z,areas,centroids,normals]=makeEllipsoidMesh(xc, ...
																  yc,zc,xr,yr,zr,Np)

[xs,ys,zs]=ellipsoid(xc,yc,zc,xr,yr,zr,Np);
[t,x,y,z]=makefliptri(xs,ys,zs); 
[tri,areas,centroids,normals]=prunemesh(t,x,y,z);
normals = -normals;  % outward pointing normals for all point based codes.