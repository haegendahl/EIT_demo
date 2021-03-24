function [J]=jacobian3dnodeFast(Node,Element,WW,Agrad,U,U0,style);

%function [J]=jacobian3d(Node,Element,WW,Agrad,U,U0,style);
% computes the Jacobian for 3D EIT when elementwise basis
% is used.
%
% Node=node data
% Element=element data
% Agrad=\int_{Element(ii) \nabla\phi_i\cdot\nabla\phi_j
% U=voltages of the injected currents
% U0=voltages of the measurement field
% style= either 'real' for reconstructing resistivity or 'comp'
% for reconstructin admittivity 

% M. Vauhkonen 31.1.2000, Univeristy of Kuopio, Finland 
% Sparse indexing modification 1.8.2006, K. Karhunen

J = 1;

%gN=max(size(Node));
%Hn=max(size(Element));

gN=max(size(Node));
HN=max(size(Element));
g=reshape([Node.Coordinate],3,gN)'; %Nodes
H=reshape([Element.Topology],4,HN)';
mH=max(max(H));
Wn=min(size(WW));

if nargin<4 
  ll = 0;
  for jj=1:gN, 
    ll = ll + length(Node(jj).ElementConnection);
  end
  Arow = zeros(ll*16,1);
  Acol = zeros(ll*16,1);
  Aval = zeros(ll*16,1);
  
  rid = 1;cid = 1;ll=0;
  for jj=1:mH
    El=Node(jj).ElementConnection;
    for ii=1:max(size(El))
      ind=Element(El(ii)).Topology; % Indices of the element
      gg=g(ind,:);
      indsig=ind;
      I=find(jj==indsig);
      anis=tedranode(gg,I);
      id = ind(:);
      idr = [id id id id];
      idc = [ind;ind;ind;ind];
      Arow(rid:rid+15,1) = idr(:) + (idc(:)-1)*gN;
      Aval(rid:rid+15,1) = anis(:);
      rid = rid + 16;      
    end     
    ll = length(El)*16;
   
    Acol(cid:cid+ll-1,1) = jj*ones(ll,1,'uint32');
    cid = cid + ll;
  end
  % Gradients of the basis functions integrated over each element. 
  J = sparse(Arow,Acol,Aval,gN^2,gN);
else
if style=='comp'
  J=zeros(size(U,2)*size(U0,2),2*size(Agrad,2));
 for ii=1:size(Agrad,2);
  Agrad1=reshape(Agrad(:,ii),gN,gN);
  JJ=-U0.'*[Agrad1,zeros(size(Agrad1));zeros(size(Agrad1)),Agrad1]*U;
  JJ=JJ(:);
  J(:,ii)=JJ;
 end,
   for ii=1+size(Agrad,2):2*size(Agrad,2);
    Agrad1=reshape(Agrad(:,ii-size(Agrad,2)),gN,gN);
    JJ=-U0.'*[zeros(size(Agrad1)),-Agrad1;Agrad1,zeros(size(Agrad1))]*U;
    JJ=JJ(:);
    J(:,ii)=JJ;
   end
elseif style=='real' 
  str = 'Patched version!';
  J=zeros(size(U,2)*size(U0,2),size(Agrad,2));
  U0 = -U0.';
  if size(U0,2)<size(U,2)
    disp([str ' loop 1'])
    for ii=1:size(Agrad,2);
      AA = reshape(Agrad(:,ii),gN,gN);
      JJ = (U0*AA)*U;      
      J(:,ii) = JJ;
    end       
  else
    disp([str ' loop 2'])
    AA = reshape(Agrad(:,1),gN,gN);
    for ii=1:size(Agrad,2);      
      JJ = U0*(reshape(Agrad(:,ii),gN,gN)*U);
      JJ = JJ(:);
      J(:,ii) = JJ;
    end
  end
end
end
