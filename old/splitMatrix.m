function [a11, a12, a21, a22] = splitMatrix(A, n1, n2);
if size(A,1) ~= n1+n2
  fprintf('splitMatrix needs dimension to be n1+n2!\n');
  a11 = 0; a12 = 0; a21 = 0; a22 = 0;
  return;
end

a11 = A(1:n1,1:n1);
a12 = A(1:n1,n1+1:end);
a21 = A(n1+1:end,1:n1);
a22 = A(n1+1:end,n1+1:end);