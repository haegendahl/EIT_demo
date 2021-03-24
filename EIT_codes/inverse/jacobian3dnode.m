function [J]=jacobian3d(Node,Element,WW,Agrad,U,U0,style);

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
Agrad=sparse(gN^2,Wn); % Gradients of the basis functions integrated over 
                       % each element. 

for jj=1:mH
Aa=sparse(gN,gN);
El=Node(jj).ElementConnection;
 for ii=1:max(size(El))
   ind=Element(El(ii)).Topology; % Indices of the element
   gg=g(ind,:);
   indsig=ind;
   I=find(jj==indsig);
   anis=tedranode(gg,I);
   Aa(ind,ind)=Aa(ind,ind)+anis;
 end
  IND = find(WW(jj,:));
  Agrad(:,IND) = Agrad(:,IND) + Aa(:);
end
J = Agrad;
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
 J=zeros(size(U,2)*size(U0,2),size(Agrad,2));
 for ii=1:size(Agrad,2);
  JJ=-U0.'*reshape(Agrad(:,ii),gN,gN)*U;
  JJ=JJ(:);
  J(:,ii)=JJ;
 end
end
end



