function [tri, x, y, z] = maketri(X,Y,Z)

N =size(X,1);
x = reshape(X, N^2,1);
y = reshape(Y, N^2,1);
z = reshape(Z, N^2,1);

tri = [];
for i=1:N-1
  for j=1:N-1
	 tri=addTriangle(tri,N,[i j; i+1 j+1; i+1 j]);
	 tri=addTriangle(tri,N,[i j; i j+1; i+1 j+1]);
  end
end
