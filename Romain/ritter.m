function foo = ritter(n, a, b, c)

% between Eqs 2.2 and 2.3.  note this goes singular as c->a
k = sqrt((a^2-b^2)/(a^2-c^2));
kprime = sqrt(1-k^2);

if ~mod(n,2)  % n even: N = n/2 types 1 (K functions), 5 (L
 % functions), 6 (M functions ), 7(N functions) in Ritter table 2

N = n/2; M = N;  % unfortunate re-use of M and N, here as counting
a1 = -(2*N -2*(1:M) +2).*(2*N+2*(1:M)-1)*k^2;
b1 = 2*N*(2*N+1)*k^2 - 4*(k^2-kprime^2)*(0:M).^2;
c1 = -(2*(0:M-1)+2).*(2*(0:M-1)+1)*kprime^2;
A1 = diag(b1)+diag(a1,-1)+diag(c1,1);
[V1,D1]= eig(A1);

M = N-1;
a5 = -(2*N-2*(1:M)).*(2*N+2*(1:M)+1)*k^2;
b5 = 2*N*(2*N+1)*k^2- (2*(0:M)+2).^2*k^2+(2*(0:M)+1).^2*kprime^2;
c5 = -(2*(0:M-1)+2).*(2*(0:M-1)+3)*kprime^2;
A5 = diag(b5)+diag(a5,-1)+diag(c5,1);
[V5,D5]= eig(A5);

a6 = -(2*N-2*(1:M)).*(2*N+2*(1:M)+1)*k^2;
b6 = 2*N*(2*N+1)*k^2- (2*(0:M)+1).^2*k^2+(2*(0:M)+2).^2*kprime^2;
c6 = -(2*(0:M-1)+2).*(2*(0:M-1)+3)*kprime^2;
A6 = diag(b6)+diag(a6,-1)+diag(c6,1);
[V6,D6]= eig(A6);

a7 = -(2*N-2*(1:M)).*(2*N+2*(1:M)+1)*k^2;
b7 = 2*N*(2*N+1)*k^2- (2*(0:M)+1).^2*(k^2-kprime^2);
c7 = -(2*(0:M-1)+2).*(2*(0:M-1)+1)*kprime^2;
A7 = diag(b7)+diag(a7,-1)+diag(c7,1);
[V7,D7]= eig(A7);

else % n odd: N = (n-1)/2
N = (n-1)/2;
a2 = -(2*N -2*(1:N) +2).*(2*N+2*(1:N)+1)*k^2;
b2 = (2*N+1)*(2*N+2)*k^2 - 4*k^2*(0:N).^2 + (2*(0:N)+1).^2*kprime^2;
c2 = -(2*(0:N-1)+2).*(2*(0:N-1)+1)*kprime^2;

A2 = diag(b2)+diag(a2,-1)+diag(c2,1);

M = N-1;
a3 = -(2*N-2*(1:M)+2).*(2*N+2*(1:M)+1)*k^2;
b3 = (2*N+1)*(2*N+2)*k^2- (2*(0:M)+1).^2*(k^2-kprime^2);
c3 = -(2*(0:M-1)+2).*(2*(0:M-1)+3)*kprime^2;n
A3 = diag(b3)+diag(a3,-1)+diag(c3,1);

a4 = -(2*N-2*(1:M)+2).*(2*N+2*(1:M)+1)*k^2;
b4 = (2*N+1)*(2*N+2)*k^2- (2*(0:M)+1).^2*k^2+4*(0:M).^2*kprime^2;
c4 = -(2*(0:M-1)+2).*(2*(0:M-1)+1)*kprime^2;
A4 = diag(b4)+diag(a4,-1)+diag(c4,1);

a8 = -(2*N-2*(1:M)).*(2*N+2*(1:M)+3)*k^2;
b8 = (2*N+1)*(2*N+2)*k^2 - (2*(0:M)+2).^2*(k^2-kprime^2);
c8 = -(2*(0:M-1)+2).*(2*(0:M-1)+3)*kprime^2;
A8 = diag(b8)+diag(a8,-1)+diag(c8,1);

end

