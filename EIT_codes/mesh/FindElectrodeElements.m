function [eltetra,E] = FindElectrodeElements(H,elnodes,Nel,order)
% [eltetra E] = FindElectrdeElements(H,elnodes,Nel,order)
%
%   elnodes{l} = indices to nodes under electrode l (cell array, size of Nel)
%   eltetra{l} = indices to elements under electrode l (cell array, size of Nel)
%   E(ii,:)    = face indices if element ii is under some electrode, zero otherwise

% A. Nissinen


% nJ=3 (1st order elements), nJ=6 (2nd order elements)
nJ = 3;  % default
if nargin==4,
  if order==2,
    nJ = 6;
  elseif order~=1
    disp('order not supported, using default')
  end
end
  

eltetra = cell(Nel,1);

nC = size(H,2);
nH = size(H,1);
J = zeros(nH,1);
E = zeros(size(H(:)),'uint32');
for el=1:Nel
  In = elnodes{el}; % nodejen indeksit elektrodilla
  J(:) = 0;
  for k=1:length(In)
    [Ir Ic] = find(H==In(k));
    J(Ir) = J(Ir)+1;
    E((Ic-1)*nH+Ir) = In(k);
  end
  eltetra{el} = uint32(find(J==nJ));
end
E = sort(reshape(E,size(H)),2);
E = E(:,nC-nJ+1:nC);
