function [Node]=MakeNode3dSmallFast(H,g);

% Function [Node]=MakeNode3dSmall(Element,g);
% computes the Node data for MeshData.
% Node is a structure including all the nodal coordinates and
% for each node there is information to which nodes (NodeConnection) 
% and elements (ElementConnection) the node is
% connected.  

% Original: M. Vauhkonen, University of Kuopio, Finland, 11.8.1999 
% Fast version with Cell-arrays by K. Karhunen, University of Kuopio, Finland, 6.7.2006

[rg,cg]=size(g);
msE = size(H,1);
rowlen = size(H,2);

maxlen = 10; % don't change
econn = zeros(rg,maxlen+1);
econn(:,1) = 1;
for k=1:msE
  id = H(k,:);
  idlen = econn(id,1);
  if max(idlen)==maxlen
    maxlen = maxlen + 10;
    swap = zeros(rg,maxlen+1);
    swap(:,1:maxlen-9) = econn;
    econn = swap;
  end
  econn(id,1) = idlen + 1;
  for ii=1:rowlen
    econn(id(ii),idlen(ii)+1) = k;
  end  
end
clear swap 

c = num2cell(g,2);
Node = cell2struct(c,'Coordinate',2);
clear c

for k=1:rg
 elen = econn(k,1); 
 Node(k).ElementConnection = [econn(k,2:elen)];
end
