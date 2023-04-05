alpha = sort(normrnd(10,2,3,1),'descend');
a = alpha(1);
b = alpha(2);
c = alpha(3);
fprintf('a=%f\tb=%f\tc=%f\n',a,b,c);
gamma = computeRomainNormalizationConstants(2, a, b, c);
gamma = gamma';

h1 = sqrt(b^2-c^2);
h2 = sqrt(a^2-c^2);
h3 = sqrt(a^2-b^2);

% Dassios eq B14
firstTerm = (a^2+b^2+c^2)/3;
secondTerm = sqrt((a^4-b^2*c^2)+(b^4-a^2*c^2)+(c^4-a^2*b^2))/3; 
LambdaD = firstTerm + secondTerm;
LambdaDprime = firstTerm - secondTerm;

analytical = [4*pi;
	      4*pi/3 *h2^2*h3^2;
	      4*pi/3 *h1^2*h3^2;
	      4*pi/3 *h1^2*h2^2;
	      -8*pi/5 * (LambdaD - LambdaDprime)*(LambdaD-a^2)* ...
	      (LambdaD-b^2)*(LambdaD-c^2);
	      8*pi/5 * (LambdaD - LambdaDprime)*(LambdaDprime-a^2)* ...
	      (LambdaDprime-b^2)*(LambdaDprime-c^2);
	      4*pi/15 *h1^2*h2^2*h3^2 *h3^2;
	      4*pi/15 *h1^2*h2^2*h3^2 *h2^2;
	      4*pi/15 *h1^2*h2^2*h3^2 *h1^2;];
	      

fprintf('Norm (analytical-numerical) = %f\n',norm(sort(analytical)- ...
						  sort(gamma)));
fprintf('Relative error = %f\n',norm(sort(analytical)- sort(gamma))/ ...
	norm(analytical));
