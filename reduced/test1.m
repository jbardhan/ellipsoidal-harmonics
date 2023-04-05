N = 20;
maxOrder = 8;
k=ceil(rand(1)*maxOrder);

V = rand(N); Vinv=inv(V);
D2=diag([logspace(3,-1,k) logspace(-1,-5,N-k)]);
A=V*D2*Vinv;


x=rand(N,1)-0.5;

Ixterm = 0*x;
for i=1:k
  Ixterm = Ixterm+V(:,i)*Vinv(i,:)*x;
end
Ixterm = Ixterm+(eye(N)-V(:,1:k)*Vinv(1:k,:))*x;

Axterm = 0*x;
for i=1:k
  Axterm = Axterm+V(:,i)*D2(i)*Vinv(i,:)*x;
end
%Axterm = Axterm+(eye(N)-V(:,1:k)*Vinv(1:k,:))*x;
Axcheck = A*x;
