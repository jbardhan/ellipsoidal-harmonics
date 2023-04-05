addpath('../Dassios');
addpath('../common');
addpath('../elliptic_package');

Ntest = [8 10 12 15];
alpha = sort((rand(1,3)*8)+1,'descend'); 

for i=1:length(Ntest)
  [tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0, alpha(1),alpha(2), ...
				      alpha(3),Ntest(i));
  filename = sprintf('../geometry/ellipsoid_points/mesh_test2_res%d.dat',Ntest(i));
  MESH = fopen(filename,'w');
  fprintf(MESH, '#  %d  %f  %f  %f\n', length(a), alpha(1),alpha(2), ...
	alpha(3));
  fprintf(MESH, '%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n', [c n a']');
  fclose(MESH);
end