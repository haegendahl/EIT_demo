function [Element,minnode] = FindElementNode(g,H,Node,d)
%function [dipind]=dipelem3(g,H,d);
% TEMPORARY VERSION NOT QUITE CORRECT
% finds element on which the dipole is.
% d has form [Dx,Dy,Dz ; rox,roy,roz]
% dipole separation is infit. small so
% only location (rox,roy,roz) is checked 
% calls function isinside which checks if
% the point is inside tetrahedra

nd = size(d,1);
Element = zeros(nd,1,'uint32');
minnode = zeros(nd,1,'uint32');

one = 1:4:16;
two = 2:4:16;
tri = 3:4:16;
four= 4:4:16;
sfl = [1:4 2:4 1 1 3:4 2 1:2 4 3];

gmax = max(g(:))*1e3;

for k=1:nd
  x = d(k,1); 
  y = d(k,2); 
  z = d(k,3);
  
  Hind = [];
  gmin = sqrt((g(:,1)-x).^2 + (g(:,2)-y).^2 + (g(:,3)-z).^2);
  [gmin gind] = min(gmin);

  %haetaan pisteessa kiinni olevat tetraedrit
  Hind = Node(gind).ElementConnection;


  %tarkistetaan onko piste elementin sisalla
  ii=0;
  yes=0;
  nHi = length(Hind);
  dface = gmax;
  while (ii<nHi) && ~yes
    ii = ii + 1;
    gid = H(Hind(ii),:);
    gid = gid(sfl);
    X = g(gid,:);
    
    %%%%%%%%%%%%%%%%%
    x1 = X(one,1);
    x2 = X(one,2);
    x3 = X(one,3);
    x21 = X(two,1)-x1;
    x31 = X(tri,1)-x1;
    y21 = X(two,2)-x2;
    y31 = X(tri,2)-x2;
    z21 = X(two,3)-x3;
    z31 = X(tri,3)-x3;

    A = y21.*z31 - y31.*z21;
    B = -(x21.*z31 - x31.*z21);
    C = x21.*y31 - x31.*y21;
    D = -x1.*A - x2.*B - x3.*C;

    sqr = sqrt(A.*A+B.*B+C.*C);
    dp4=(A.*X(four,1)+B.*X(four,2)+C.*X(four,3)+D)./sqr;
    dp=(A*x+B*y+C*z+D)./sqr;

    check = dp.*dp4 >= 0;
    if all(check), yes = 1; end
    %%%%%%%%%%    
  end
  if yes, Element(k,1) = Hind(ii); else, minnode(k,1) = gind; end
end

