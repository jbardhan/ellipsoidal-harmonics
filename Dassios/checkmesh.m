function [areas,centroids,normals,bad,hand]= checkmesh(tri,x,y,z)

points = [x y z];
hand = 0;
for i=1:size(tri,1)
  V = points(tri(i,:),:)';
  e1 = V(:,2) - V(:,1); 
  e2 = V(:,3) - V(:,1); 
  areas(i) = 0.5*norm(cross(e1,e2));
  e1=e1/norm(e1);
  e2=e2/norm(e2);
  normals(:,i) = cross(e1,e2);
  normals(:,i) = normals(:,i)/norm(normals(:,i));
  centroids(:,i) = mean(V,2);
  hand(i)=normals(:,i)'*centroids(:,i);
  if norm(hand(i)) < 1e-10
	 bad(i) = 1;
  else
	 bad(i) = 0;
  end
end

normals = normals';
centroids = centroids';
