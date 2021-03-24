function [H,Element,nc] = organize_tetras2nd(H,Element,g)
% Netgen & Matlab (own meshing to 2nd degree) safe

nE = size(H,1);

nc = 0;
for ii=1:nE
  nodes = H(ii,:);  % tetra nodes
  mp = nodes(1:4);  % main nodes (this may fail with exotic mesh formats)
  cp = .5*(g(mp([1 2 1 1 2 3]),:) + g(mp([2 3 3 4 4 4]),:));
  gg = g(nodes(5:10),:);  % center nodes (2nd order)
  for jj=1:6
    [m I(jj)] = min((cp(jj,1)-gg(:,1)).^2 + (cp(jj,2)-gg(:,2)).^2 + (cp(jj,3)-gg(:,3)).^2);
  end
    
  nodes2 = [mp nodes(I+4)];
  H(ii,:) = nodes2;
  Element(ii).Topology = nodes2;
  if max(nodes - nodes2), nc = nc + 1; end  
end
