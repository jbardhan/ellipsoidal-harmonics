addpath('../Dassios');
addpath('../common');
addpath('../elliptic_package');

Ntest = [5 8 10 12 15];
alpha = (rand(1,3)*8)+1; 
numtosave = 10;
f=@(x,y) 1./((x+y^2).*sqrt((x+alpha(1)^2).*(x+alpha(2)^2).* (x+alpha(3)^2)));
analytical = (prod(alpha)*[quad(@(x)f(x,alpha(1)),0,1e6) quad(@(x)f(x, ...
						  alpha(2)),0,1e6) ...
		    quad(@(x)f(x,alpha(3)),0,1e6)]-1)/2;
for i=1:length(Ntest)
  [tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0, alpha(1),alpha(2), ...
				      alpha(3),Ntest(i));
  areas(i) = sum(a);
  Npanels(i) = length(a);

  % generate centroid based approximation to K
  [V,K]=genPointMatrices(c,n,a);
  [K_eigenvecs,K_eigenvals] = eig(K);
  kd=real(diag(K_eigenvals));
  [junk,m1]=max(abs(c(:,1)'*K_eigenvecs)); 
  [junk,m2]=max(abs(c(:,2)'*K_eigenvecs)); 
  [junk,m3]=max(abs(c(:,3)'*K_eigenvecs)); 
  Ksave(i,:) = kd([m1 m2 m3])';

  % generate BEM-based approximation to K
  meshX = x(tri)'; meshY = y(tri)'; meshZ = z(tri)';
  meshdata = struct('face',tri,'X',meshX,'Y',meshY,'Z',meshZ);
  [V,K_bem]=colloc_Laplace(meshdata,c,n,a);  

  K_bem = K_bem-diag(diag(K_bem));
  K_bem = K_bem;															
  [K_bem_eigenvecs,K_bem_eigenvals]=eig(K_bem);
  kdbem=real(diag(K_bem_eigenvals));
  [junk,m1]=max(abs(c(:,1)'*K_bem_eigenvecs)); 
  [junk,m2]=max(abs(c(:,2)'*K_bem_eigenvecs)); 
  [junk,m3]=max(abs(c(:,3)'*K_bem_eigenvecs)); 
  Kbemsave(i,:) = kdbem([m1 m2 m3])';
  
  resultsBEM(i) = struct('eigenvalues',sort(kdbem,'ascend'));
end

x = sort(alpha);
a = x(1); b = x(2); c = x(3);

nmax = 12;

lambda = a;
Imn = computeExteriorIntegral(lambda,nmax,a,b,c); % computes all

index = 1;
for i=0:nmax
  for j=0:2*i
    [Emn(index),Emnderiv(index)] = calcLame(lambda,i,j,a,b,c);
    index = index+1;
  end
end

chainRuleRitter = 1/a;
lambda = ((2*a*b*c)*(Emnderiv.*chainRuleRitter).*(Imn.*Emn) -1)/2;
