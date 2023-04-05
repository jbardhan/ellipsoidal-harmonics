for loopvar = 1:1
  vars = rand(3,1) * 10; vars=sort(vars,'descend');
  a = 7.3543; %vars(1);
  b = 5.2283; %vars(2);
  c = 3.4897; %vars(3);
  
nmax = 5;

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

index = 1;

% upper limit needs to be checked more carefully in the future.
% Python and Matlab currently have some differences because Matlab
% offers an easy adaptive quadrature and NumPy does not (yet).
intVal(1)  = quad(@(t) 1./((t+a^2).*sqrt((t+a^2).*(t+b^2).*(t+c.^2))), 0, 1e6);
dipAnal(1) = (a*b*c * intVal(1) - 1)/2;
intVal(2)  = quad(@(t) 1./((t+b^2).*sqrt((t+a^2).*(t+b^2).*(t+c.^2))), 0, 1e6);
dipAnal(2) = (a*b*c * intVal(2) - 1)/2;
intVal(3)  = quad(@(t) 1./((t+c^2).*sqrt((t+a^2).*(t+b^2).*(t+c.^2))), 0, 1e6);
dipAnal(3) = (a*b*c * intVal(3) - 1)/2;

fprintf('a = %f, b = %f, c = %f\n',a,b,c);
fprintf('analytical:'); 
fprintf('%3.5f   ',dipAnal); 
fprintf('\n');

fprintf('numerical: ');
fprintf('%3.5f   ',lambda(2:4)); 
fprintf('\n\n');

end

