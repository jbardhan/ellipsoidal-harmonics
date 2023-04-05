spherecenter= [ 0 0 0];
radius = 5;
np = [576];% 750 1000 1500 2000];
for i=1:length(np)
  [pointlist,weights,ang,normals] = getSphPoints(spherecenter, radius, np(i));
  filename = sprintf('../geometry/ellipsoid_points/sphere_mesh_rad%d_%d.dat',radius,np(i));
  a=fopen(filename,'w');
  fprintf(a,'# %d  %f %f %f\n',np(i),radius,radius,radius);
  fprintf(a, '%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n', [pointlist normals weights]');
  fclose(a);
end
