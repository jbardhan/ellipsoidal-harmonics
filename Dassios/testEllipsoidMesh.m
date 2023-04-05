addpath('../common');
addpath('../elliptic_package');

Ntest = [3 5 8 10 12 15];
test1 = 0;
test2 = 0;
test3 = 0;
test4 = 1;

if test1
% test 1: sphere area
for i=1:length(Ntest)
  r=1;
  [tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0,r,r,r,Ntest(i));
  areas(i) = sum(a);
  Npanels(i) = length(a);
end
analyticalSph = 4*pi*r^2;
figure;
a=loglog(Npanels,abs(analyticalSph-areas),'-o','linewidth',2);
set(gca,'fontsize',16);
xlabel('Number of triangles');
ylabel('Deviation from analytical');
title('Test 1: Convergence of sphere surface area');
end

if test2
% test 2: ellipsoid areas
alpha = (rand(1,3)*8)+1;
for i=1:length(Ntest)
  [tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0,alpha(1),alpha(2), ...
												  alpha(3),Ntest(i));
  areas(i) = sum(a);
  Npanels(i) = length(a);
end
alphasort=sort(alpha,'descend');
aa=alphasort(1); bb=alphasort(2); cc=alphasort(3);
alphaArea = acos(cc/aa);
m = (bb^2-cc^2)/bb^2/sin(alphaArea)/sin(alphaArea);
[F,E]= elliptic12(alphaArea,m,1e-6);
analytical=2*pi*(cc^2+bb*sqrt(aa^2-cc^2)*E + (bb*cc^2/sqrt(aa^2- ...
																  cc^2))*F);
figure;
a=loglog(Npanels,abs(analytical-areas),'-o','linewidth',2);
set(gca,'fontsize',16);
xlabel('Number of triangles');
ylabel('Deviation from analytical');
title('Test 2: Convergence of ellipsoid surface area');
end

if test3
numtosave= 9;
% test 3: eigendecomposition of K for sphere
  r=2;
  [lambdaV,lambdaK]=computeLaplaceEigenvalues(r,100);

for i=1:length(Ntest)
  [tri,x,y,z,a,c,n]=makeEllipsoidMesh(0,0,0,r,r,r,Ntest(i));
  areas(i) = sum(a);
  Npanels(i) = length(a);

  % generate centroid based approximation to K
  [V,K]=genPointMatrices(c,n,a);
  [K_eigenvecs,K_eigenvals] = eig(K);
  K_eigenvals = sort(diag(real(K_eigenvals)),'ascend');
  Ksave(i,:)=K_eigenvals(1:numtosave);

  % generate BEM-based approximation to K
  meshX = x(tri)'; meshY = y(tri)'; meshZ = z(tri)';
  meshdata = struct('face',tri,'X',meshX,'Y',meshY,'Z',meshZ);
  [V,K_bem]=colloc_Laplace(meshdata,c,n,a);
  keyboard
  K_bem = -K_bem; % calcp uses other handedness. 
  K_bem = K_bem-diag(diag(K_bem));  % zero out diagonal.  diagonal
                                    % term for flat panels = 0.
  [K_bem_eigenvecs,K_bem_eigenvals]=eig(K_bem);
  K_bem_eigenvals = sort(diag(real(K_bem_eigenvals)),'ascend');
  Kbemsave(i,:)=K_bem_eigenvals(1:numtosave);
end  
figure;
loglog(Npanels,abs(lambdaK(1:numtosave)*ones(1,length(Npanels))- ...
						 Ksave'),'--o');
hold on;
loglog(Npanels,abs(lambdaK(1:numtosave)*ones(1,length(Npanels))- ...
						 Kbemsave'),'-s');
set(gca,'fontsize',16);
xlabel('Number of triangles');
ylabel('Deviation from analytical eigenvalue');
title('Test 3: Convergence of Operator Eigenvalues');
end

if test4
% test 4: eigendecomposition of K for ellipsoid
alpha = (rand(1,3)*8)+1;
numtosave = 10;
f=@(x,y) 1./((x+y^2).*sqrt((x+alpha(1)^2).*(x+alpha(2)^2).* ...
											  (x+alpha(3)^2)));
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
  K_bem = -K_bem;															
  [K_bem_eigenvecs,K_bem_eigenvals]=eig(K_bem);
  kdbem=real(diag(K_bem_eigenvals));
  [junk,m1]=max(abs(c(:,1)'*K_bem_eigenvecs)); 
  [junk,m2]=max(abs(c(:,2)'*K_bem_eigenvecs)); 
  [junk,m3]=max(abs(c(:,3)'*K_bem_eigenvecs)); 
  Kbemsave(i,:) = kdbem([m1 m2 m3])';
  
end

figure;set(gca,'fontsize',16);
loglog(Npanels,abs(ones(length(Npanels),1)*analytical-Ksave),'--o');
hold on;
loglog(Npanels,abs(ones(length(Npanels),1)*analytical-Kbemsave),'-s');
title('Test 4: Convergence of Ellipsoid-Operator Eigenvalues');
xlabel('Number of triangles');
ylabel('Deviation from analytical');
end