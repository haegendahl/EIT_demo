function int=tetraLinSigma(g,s,ip,L);

% Lasketaan sigman arvot integrointipisteissä
sigma = [1-ip(:,1)-ip(:,2)-ip(:,3),ip(:,1),ip(:,2),ip(:,3)]*s;

Jt = L*g;
iJt = inv(Jt);
dJt = abs(det(Jt));
G = iJt*L;
GdJt = G'*G*dJt/24;
int = sum(sigma)*GdJt;
