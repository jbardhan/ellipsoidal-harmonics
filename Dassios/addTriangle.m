function tri = addTriangle(tri, N, indices)
indices(:,2)=indices(:,2)-1;
newTri = (indices * [1; N])';
tri = [tri; newTri];