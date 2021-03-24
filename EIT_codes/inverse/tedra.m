function int=tedra(g);

int = 1;

% P.Ronkanen 28.6.1996

% g=globaalin elementin solmupisteet


w=[1/24*ones(4,1)];
a=0.58541020;
b=0.13819660;
ip=[b b b;a b b;b a b;b b a];
L=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
Jt=L*g;
iJt=inv(Jt);
dJt=abs(det(Jt));
G=iJt*L;
int=1/6*G'*G*dJt;















