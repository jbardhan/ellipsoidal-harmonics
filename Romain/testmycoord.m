addpath('../Dassios');
Ntest = 10;
alpha = sort((rand(1,3)*8)+1,'descend');
[tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0,alpha(1),alpha(2), ...
												  alpha(3),Ntest);
