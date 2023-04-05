a = 3; 
b = 2;
c = 1;
h = sqrt(a^2-b^2);
k = sqrt(a^2-c^2);
%lambdatest = [8.4853 2.3872 2.2338];
lambdatest = [2*k (h+k)/2 h/2;
				  2*k -(h+k)/2 h/2;
				  2*k (h+k)/2 -h/2;
				  2*k -(h+k)/2 -h/2];

for i=1:size(lambdatest,1)
    pnt(i,:) = ellToCart(a,b,c,lambdatest(i,:))
end
