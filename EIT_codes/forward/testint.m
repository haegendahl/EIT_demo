function int = testint(g)

dx = g(2:4,1) - g(1,1);
dy = g(2:4,2) - g(1,2);
dz = g(2:4,3) - g(1,3);

Jt = [dx dy dz];
dJt = det(Jt);
adJt = abs(dJt)/(6*dJt);

iJt = [dy(2)*dz(3)-dz(2)*dy(3),dz(1)*dy(3)-dy(1)*dz(3),dy(1)*dz(2)-dy(2)*dz(1);
       dz(2)*dx(3)-dx(2)*dz(3),dx(1)*dz(3)-dx(3)*dz(1),dz(1)*dx(2)-dx(1)*dz(2);
       dx(2)*dy(3)-dx(3)*dy(2),dx(3)*dy(1)-dx(1)*dy(3),dx(1)*dy(2)-dy(1)*dx(2)];
G = [-sum(iJt(1,:)) iJt(1,:);
     -sum(iJt(2,:)) iJt(2,:);
     -sum(iJt(3,:)) iJt(3,:)]*adJt;
int = G'*G;
