function [newtri,areas,centroids,normals] = prunemesh(t,X,Y,Z)
thresholdArea = 1e-5;

[areas,centroids,normals]=checkmesh(t,X,Y,Z);
good= find(areas>thresholdArea);

newtri = t(good,:);
areas = areas(good);
centroids = centroids(good,:);
normals = normals(good,:);

if 0 % only use if you're generating ellipsoid mesh with the ugly
     % "fromEllipsoidalCoord" function... which you shouldn't.
for i=1:length(good)
  notduplicated(i) = 1;
  for j=1:i-1
	 distCent = norm(centroids(i,:)-centroids(j,:));
	 if distCent < 1e-8
		notduplicated(j) = 0;
	 end
  end
end

good = find(notduplicated);
newtri = t(good,:);
areas = areas(good);
centroids = centroids(good,:);
normals = normals(good,:);
end