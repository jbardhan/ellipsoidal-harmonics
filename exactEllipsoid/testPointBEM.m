addpath('../Dassios');
addpath('../common');
addpath('../elliptic_package');

pqrData = readpqr('../geometry/ellipsoid_points/monopole.pqr');
epsIn = 4;
epsOut = 80;

nplist = [100 250 400];
for i=1:length(nplist)
  filename = sprintf('../geometry/ellipsoid_points/sphere_mesh_rad3_%d.dat',nplist(i));

  FID = fopen(filename, 'r');
  numpanels = fscanf(FID,'# %d', 1);
  headerdata = fscanf(FID,'%f ', 3);
  a = headerdata(1); b = headerdata(2); c = headerdata(3);
  frewind(FID);
  pointdataCell = textscan(FID,'%f %f %f %f %f %f %f','Headerlines',1);
  pointdata = cell2mat(pointdataCell);
  fclose(FID);
  points = pointdata(:,1:3);  % really centroids
  normals = pointdata(:,4:6);
  areas = pointdata(:,7);
  
  [A,B,C] = getPointECF(pqrData, points, normals, areas, epsIn, ...
				 epsOut);
keyboard
  energylist(i) = 0.5 * 332.112 * C * (A\B);
end

analytical = 0.5 * 332.112 * (1/epsOut - 1/epsIn) * 1/a;

figure;
set(gca,'fontsize',16);
H = loglog(nplist, abs(energylist-analytical),'s--','linewidth',2,'markersize',12);
xlabel('Number of points');
ylabel('Absolute error (kcal/mol)');



