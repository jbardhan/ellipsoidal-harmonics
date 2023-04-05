addpath('../Dassios');
addpath('../common');
addpath('../elliptic_package');

pqrData = readpqr('../geometry/ellipsoid_points/monopole.pqr');
epsIn = 4;
epsOut = 80;

reslist = [15];
for i=1:length(reslist)

  % step 1. load data
  filename = sprintf('../geometry/ellipsoid_points/mesh_test1_res%d.dat',reslist(i));
  FID = fopen(filename, 'r');
  numpanels = fscanf(FID,'# %d', 1);
  headerdata = fscanf(FID,'%f ', 3);
  a = headerdata(1); b = headerdata(2); c = headerdata(3);
  frewind(FID);
  pointdataCell = textscan(FID,'%f %f %f %f %f %f %f','Headerlines',1);
  pointdata = cell2mat(pointdataCell);
  fclose(FID);

  % step 2. generate relevant geometry variables and operators
  points = pointdata(:,1:3);  % really centroids of a mesh,
                              % i.e. not exactly on the ellipsoid surface.
  normals = pointdata(:,4:6);  
  areas = pointdata(:,7);
  [A,B,C,K] = getPointECF(pqrData, points, normals, areas, epsIn, ...
				 epsOut);

  % step 3. compute energy 
  energylist(i) = 0.5 * 332.112 * pqrData.q' * C * (A\(B*pqrData.q));
  numpanellist(i) = length(areas);
end

