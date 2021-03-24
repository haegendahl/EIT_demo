function int=tedra(g,I);

% P.Ronkanen 28.6.1996

% g=globaalin elementin solmupisteet
int = 1;


w=[1/24*ones(4,1)];
a=0.58541020;
b=0.13819660;
ip=[b b b;a b b;b a b;b b a];
L=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
Jt=L*g;
iJt=inv(Jt);
dJt=abs(det(Jt));
G=iJt*L;
GdJt=G'*G*dJt;

S=[[1-ip(1,1)-ip(1,2)-ip(1,3);ip(1,1);ip(1,2);ip(1,3)],[1-ip(2,1)-ip(2,2)-ip(2,3);ip(2,1);ip(2,2);ip(2,3)], ...
   [1-ip(3,1)-ip(3,2)-ip(3,3);ip(3,1);ip(3,2);ip(3,3)],[1-ip(4,1)-ip(4,2)-ip(4,3);ip(4,1);ip(4,2);ip(4,3)]]*w;
 int=S(I)*GdJt;
















