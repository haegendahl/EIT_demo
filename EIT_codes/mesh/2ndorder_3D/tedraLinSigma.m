function int=tedraLinSigma(g,sigma);
%
w=[1/24*ones(4,1)];
a=0.58541020;
b=0.13819660;
ip=[b b b;a b b;b a b;b b a];
%g=g([3 1 7 5 2 8 9 4 6 10],:);
int=0;
for ii=1:4
k=ip(ii,1);
n=ip(ii,2);
c=ip(ii,3);
L=[4*k-3+4*n+4*c 4*k-1 0 0 -8*k+4-4*n-4*c 4*n -4*n -4*c 4*c 0;
   4*n-3+4*k+4*c 0 4*n-1 0 -4*k 4*k -8*n+4-4*k-4*c -4*c 0 4*c;
   4*c-3+4*k+4*n 0 0 4*c-1 -4*k 0 -4*n -8*c+4-4*k-4*n 4*k 4*n];
S=[1-k-n-c;k;n;c];
 Jt=L*g;
 iJt=inv(Jt);
 dJt=abs(det(Jt));
 G=iJt*L;
 int=int+w(ii)*S'*sigma*dJt*G'*G;
end
