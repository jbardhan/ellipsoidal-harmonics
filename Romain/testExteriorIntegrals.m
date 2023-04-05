a = 3;
b = 2;
c = 1;

h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);

lambda = a*2;
I = computeExteriorIntegral(lambda, 2,a,b,c);

rootSum = sqrt((a^4-b^2*c^2)+(b^4-a^2*c^2)+(c^4-a^2*b^2));
DassiosLambda =      (1/3)*(a^2+b^2+c^2) + (1/3)*rootSum;
DassiosLambdaPrime = (1/3)*(a^2+b^2+c^2) - (1/3)*rootSum;

fprintf('Dassios D1 error:\n');
lhs = 3 * (DassiosLambda+DassiosLambdaPrime);
rhs = 2 * (a^2 + b^2 + c^2);
fprintf('%f = %f - %f\n\n', lhs-rhs, lhs, rhs);

fprintf('Dassios D2 error:\n');
lhs = 3 * DassiosLambda * DassiosLambdaPrime;
rhs = a^2*b^2 + a^2*c^2 + b^2*c^2;
fprintf('%f = %f - %f\n\n', lhs-rhs, lhs, rhs);

fprintf('Dassios D3 error:\n');
lhs = ((-1)*(b^2-c^2)*(DassiosLambda-a^2)) + ...
      ((-1)^2*k^2*(DassiosLambda-b^2)) + ...
      ((-1)^3*h^2*(DassiosLambda-c^2));
rhs =  ((-1)*(b^2-c^2)*(DassiosLambdaPrime-a^2)) + ...
      ((-1)^2*k^2*(DassiosLambdaPrime-b^2)) + ...
      ((-1)^3*h^2*(DassiosLambdaPrime-c^2));
fprintf('%f should = 0\n',lhs);
fprintf('%f should = 0\n\n',rhs);

fprintf('Dassios D4 error:\n');
lhs = ((-1)^1*a^2*(b^2-c^2)*(DassiosLambda-a^2)) + ...
      ((-1)^2*b^2*k^2*(DassiosLambda-b^2)) + ...
      ((-1)^3*c^2*h^2*(DassiosLambda-c^2));
rhs = ((-1)^1*a^2*(b^2-c^2)*(DassiosLambdaPrime-a^2)) + ...
      ((-1)^2*b^2*k^2*(DassiosLambdaPrime-b^2)) + ...
      ((-1)^3*c^2*h^2*(DassiosLambdaPrime-c^2));
analytical = h^2*k^2*(b^2-c^2);
fprintf('%f = %f (numerical) - %f (analytical)\n',lhs-analytical, ...
	lhs,analytical);
fprintf('%f = %f (numerical) - %f (analytical)\n\n',rhs-analytical, ...
	rhs,analytical);

fprintf('Dassios D5 error:\n');
alpha = [a b c];
lhs = sum(alpha.^2 ./(alpha.^2 - DassiosLambda));
rhs = sum(alpha.^2 ./(alpha.^2 - DassiosLambdaPrime));
fprintf('%f = %f (numerical) - 3 (analytical)\n',lhs-3, lhs);
fprintf('%f = %f (numerical) - 3 (analytical)\n\n',rhs-3, rhs);
clear alpha;

fprintf('Dassios D6 error:\n');
anal1 = (-1)^(1+1)*h^2*k^2*(b^2-c^2);
num1  = 3*(b^2-c^2)*(DassiosLambda - a^2)*(DassiosLambdaPrime-a^2);

anal2 = (-1)^(2+1)*h^2*k^2*(b^2-c^2);
num2  = 3*k^2*(DassiosLambda - b^2)*(DassiosLambdaPrime-b^2);

anal3 = (-1)^(3+1)*h^2*k^2*(b^2-c^2);
num3  = 3*h^2*(DassiosLambda - c^2)*(DassiosLambdaPrime-c^2);
fprintf('%f = %f (numerical) - %f (analytical)\n',num1-anal1,num1, ...
	anal1);
fprintf('%f = %f (numerical) - %f (analytical)\n',num2-anal2,num2, ...
	anal2);
fprintf('%f = %f (numerical) - %f (analytical)\n\n',num3-anal3,num3, ...
	anal3);


% Dassios D7: their h_3 = our h and their h_2 = our k 
analyticalSumTest1 = 1.0/(lambda*sqrt(lambda^2-h^2)*sqrt(lambda^2-k^2));
numericalSumTest1 = sum(I(2:4));

fprintf('Dassios D7 error:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest1-analyticalSumTest1, numericalSumTest1, ...
	analyticalSumTest1);

% Dassios D8: their alpha1 = our a, their alpha2 = our b, their alpha3
% ... = our c
analyticalSumTest2 = I(1) - (lambda^2-a^2)/(lambda*sqrt(lambda^2- ...
						  h^2)*sqrt(lambda^2-k^2));
numericalSumTest2 = a^2*I(2) + b^2*I(3) + c^2 * I(4);

fprintf('Dassios D8 error:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest2-analyticalSumTest2,numericalSumTest2, ...
	analyticalSumTest2);

% Dassios D9+D10: their Lambda = our DassiosLambda, their LambdaPrime
% ... = our DassiosLambdaPrime
analyticalSumTest3 = 1/(2*(DassiosLambda-a^2+lambda^2)*lambda* ...
			sqrt(lambda^2-h^2)*sqrt(lambda^2-k^2)) - ...
    (1.0/2.0)*(I(2)/(DassiosLambda-a^2) + I(3)/(DassiosLambda-b^2) ...
	       + I(4)/(DassiosLambda-c^2));
numericalSumTest3 = I(5);

fprintf(['The following two tests are broken.\nCorrectness depends' ...
	 ' on a different normalization of the Lame functions.\n\n']); 
fprintf('Dassios D9 error is broken:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest3-analyticalSumTest3,numericalSumTest3, ...
	analyticalSumTest3);

analyticalSumTest4 = 1/(2*(DassiosLambdaPrime-a^2+lambda^2)*lambda* ...
			sqrt(lambda^2-h^2)*sqrt(lambda^2-k^2)) - ...
    (1.0/2.0)*(I(2)/(DassiosLambdaPrime-a^2) + I(3)/ ...
	       (DassiosLambdaPrime-b^2) + I(4)/(DassiosLambdaPrime-c^2));
numericalSumTest4 = I(6);

fprintf('Dassios D10 error is broken:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest4-analyticalSumTest4,numericalSumTest4, ...
	analyticalSumTest4);

% Dassios D11,D12,D13: their h1^2 = our (b^2-c^2)
analyticalSumTest5 = (1.0/h^2) * (I(3)-I(2));
numericalSumTest5  = I(7);

fprintf('Dassios D11 error:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest5-analyticalSumTest5,numericalSumTest5, ...
	analyticalSumTest5);

analyticalSumTest6 = (1.0/k^2) * (I(4)-I(2));
numericalSumTest6  = I(8);

fprintf('Dassios D12 error:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest6-analyticalSumTest6,numericalSumTest6, ...
	analyticalSumTest6);

analyticalSumTest7 = (1.0/(b^2-c^2)) * (I(4)-I(3));
numericalSumTest7  = I(9);

fprintf('Dassios D13 error:\n');
fprintf('%f = %f (numerical) - %f (analytical)\n\n', ...
	numericalSumTest7-analyticalSumTest7,numericalSumTest7, ...
	analyticalSumTest7);

