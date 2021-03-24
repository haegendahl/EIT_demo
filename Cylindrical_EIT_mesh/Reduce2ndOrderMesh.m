function [H,g,Node,Element,elind,eltetra,E] = Reduce2ndOrderMesh(H2,g2,elind2,meshname)
% Reduces mesh from 2nd to 1st order basis
% Netgen compatible mesh format expected (linear elements are at the top)

Nel = length(elind2);

H = H2(:,1:4);
ng = max(H(:));
g = g2(1:ng,:);

elind = cell(Nel,1);
for ii=1:Nel
  I = elind2{ii};
  J = find(I<=ng);
  elind{ii} = I(J);
end

Node = MakeNode3dSmallFast(H,g);
[eltetra,E] = FindElectrodeElements2(H,Node,elind,1);
Element = MakeElement3dSmallCellFast(H,eltetra,E);

eval(['save ',meshname,' g H Node Element eltetra E'])







