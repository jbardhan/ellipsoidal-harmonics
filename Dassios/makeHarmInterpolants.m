function V = makeHarmInterpolants(maxorder,grid,X,Y,Z)

N =size(X,1);
x = reshape(X, N^2,1);
y = reshape(Y, N^2,1);
z = reshape(Z, N^2,1);
MU = reshape(grid.MU,N^2,1);
NU = reshape(grid.NU,N^2,1);

count = 1;
V = zeros(count,9); % maxorder = 2 --> 9 
tri = [];
for i=1:N-1
  for j=1:N-1
	 triVerts=addTriangle(tri,N,[i j; i+1 j+1; i+1 j]);
	 V(count,:) = evalHarmonics(maxorder,grid,MU,NU,triVerts);
	 count = count+1;
	 triVerts=addTriangle(tri,N,[i j; i j+1; i+1 j+1]);
	 V(count,:) = evalHarmonics(maxorder,grid,MU,NU,triVerts);
	 count=count+1;
  end
end
