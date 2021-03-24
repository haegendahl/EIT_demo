function [eltetra,E,nc] = FindElectrodeElements2(H,Node,elnodes,order,gindx)
% [eltetra E] = FindElectrodeElements2(H,Node,elnodes,order,gindx)
%
%   elnodes{l} = indices to nodes under electrode l (cell array, size of Nel)
%   eltetra{l} = indices to elements under electrode l (cell array, size of Nel)
%   E(ii,:)    = face indices if element ii is under some electrode, zero otherwise
%   gindx      = (optional) reindexing basing on geometry, 'yes' or 'no'. 
%                  This has no effect if 1st order basis is used. Default = 'yes'
%   nc         = how many indices were changed (related to gindx)


nc = 0;

% nJ=3 (1st order elements), nJ=6 (2nd order elements)
nJ = 3;  % default
if nargin>3,
  if order==2,
    nJ = 6;
  elseif order~=1
    disp('order not supported, using default')
  end
end
if nargin<5, gindx = 'yes'; end

nH = size(H,1);

Nel = prod(size(elnodes));
if size(elnodes,2)>1
  elnodes = reshape(elnodes,Nel,1);
end


E = zeros(nH,nJ,'uint32');
eltetra = cell(Nel,1);

tetra_mask = zeros(nH,1,'uint8');

% loop through every electrode
for ii=1:Nel
  ind_node = elnodes{ii};  
  node_len = length(ind_node);
  
  % loop through nodes in electrode ii
  for jj=1:node_len
    ind_tetra = Node(ind_node(jj)).ElementConnection;    
    tetra_len = length(ind_tetra);
    
    % check every tetrahedron connected to node ptr *ii
    for kk=1:tetra_len
      tetra = ind_tetra(kk);
      Hind = H(tetra,:);
      
      if ~tetra_mask(tetra)  % ...we want select the tetrahedron only once
        [C,II] = intersect(Hind,ind_node);
        if (length(C)==nJ)
          eltetra{ii} = [eltetra{ii};uint32(tetra)];
          E(tetra,:) = Hind(sort(II));
        end      
      end
      tetra_mask(tetra)=1;
    end
    
  end
  
end

if (strcmpi(gindx,'yes')) && (order==2), [E,nc] = reindex(E,Node); end
return

%% reindex faces basing on geometry 
%%   note: this may fail with severely skewed elements
function [E,nc] = reindex(E,Node)

gN = max(size(Node));
g = reshape([Node.Coordinate],3,gN)'; %Nodes
nE = size(E,1);

nc = 0;
for ii=1:nE
  if all(E(ii,:))
    nodes = E(ii,:);
    mp = nodes(1:3);
    cp = .5*(g(mp,:) + g(mp([2 3 1]),:));
    gg = g(nodes(4:6),:);  % center nodes (2nd order)
    [m I1] = min((cp(1,1)-gg(:,1)).^2 + (cp(1,2)-gg(:,2)).^2 + (cp(1,3)-gg(:,3)).^2);
    [m I2] = min((cp(2,1)-gg(:,1)).^2 + (cp(2,2)-gg(:,2)).^2 + (cp(2,3)-gg(:,3)).^2);
    [m I3] = min((cp(3,1)-gg(:,1)).^2 + (cp(3,2)-gg(:,2)).^2 + (cp(3,3)-gg(:,3)).^2);
    nodes2 = [mp nodes([I1 I2 I3]+3)];
    E(ii,:) = nodes2;
    if ~isequal(nodes,nodes2), nc = nc + 1; end
  end
end

disp(['Reindexing: ' num2str(nc) ' face(s) changed!'])
