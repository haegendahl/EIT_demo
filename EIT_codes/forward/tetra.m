function int=tetra(g,L);

Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L;
int = G'*G*dJt/24;

