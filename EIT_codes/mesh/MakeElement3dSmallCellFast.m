function [Element]=MakeElement3dSmallCellFast(H,eltetra,E);

% Function [Element]=MakeElement3dsmall(H,bg,E); 
% computes the Element data for MeshData. 
% Element is a structure including the topology (H) of each element,
% faces of each element, adjacent element of each face and information
% whether the face is internal or boundary face.
% elind includes the boundary node indices that are under the electrodes
% E includes the eleement indices that are under the electrodes

% M. Vauhkonen 10.10.1998. Modified for EIDORS3D 28.1.2000 by 
% M. Vauhkonen, University of Kuopio, Finland.

% Modified to use cell arrays 6.7.2006 by K. Karhunen, University
% of Kuopio, Finland. This modification allows arbitrary number of
% elements under each electrode.

% Speed improvement, 26.07.2006 K. Karhunen

[rH,cH]=size(H);
nel = size(eltetra,1);

c = num2cell(H,2);
Element = cell2struct(c,'Topology',2);
clear c

[Element(1:rH).Electrode] = deal([]);

for k=1:nel
  id = eltetra{k};
  for n=1:length(id)
    ii = id(n);
    Element(ii).Electrode{1}=k;
    Element(ii).Electrode{2}=[E(ii,:)];
  end  
end
 
