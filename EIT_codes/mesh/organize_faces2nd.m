function [Element2,nc] = organize_faces2nd(Element,g)
% Netgen & Matlab (own meshing to 2nd degree) safe

nE = size(Element,1);
Element2 = Element;

nc = 0;
for ii=1:nE
  if ~isempty(Element(ii).Electrode)
    nodes = Element(ii).Electrode{2};  % face nodes
    mp = nodes(1:3);  % main nodes (this may fail with exotic mesh formats)
    cp = .5*(g(mp,:) + g(mp([2 3 1]),:));
    gg = g(nodes(4:6),:);  % center nodes (2nd order)
    [m I1] = min((cp(:,1)-gg(1,1)).^2 + (cp(:,2)-gg(1,2)).^2 + (cp(:,3)-gg(1,3)).^2);
    [m I2] = min((cp(:,1)-gg(2,1)).^2 + (cp(:,2)-gg(2,2)).^2 + (cp(:,3)-gg(2,3)).^2);
    [m I3] = min((cp(:,1)-gg(3,1)).^2 + (cp(:,2)-gg(3,2)).^2 + (cp(:,3)-gg(3,3)).^2);
    nodes2 = [mp nodes([I1 I2 I3]+3)];
    Element2(ii).Electrode{2} = nodes2;
    if max(nodes - nodes2), nc = nc + 1; end
  end
end
