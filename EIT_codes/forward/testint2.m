function int = testint2(g,L,Jt)

dx = g(2:4,1) - g(1,1);
dy = g(2:4,2) - g(1,2);
dz = g(2:4,3) - g(1,3);
Jt = [dx dy dz];

iJt = inv(Jt);
dJt = 1;%abs(det(Jt))/6;
G = iJt*L;
int = G'*G*dJt;
