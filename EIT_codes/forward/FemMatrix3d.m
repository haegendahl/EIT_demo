function [A]=FemMatrix3d(Node,Element,sigma,z,C);

% function [A]=FemMatrix(Node,Element,sigma,z,C) supplies matrices
% needed for the 3D FEM solution in the EIT forward problem. 
% gn contains the nodes in (x,y) coordinates, 
% Hn the indexes to gn of the tetrahedron vertices.
% Eind contains the indices of the tetrahedrons that have
% two nodes on the boundary and are under the electrodes. 
% elind contains on each row the nodes of the elements
% that are under the coresponding electrode.
% rho is the resistivity vector and z the contact impedance
% vector.

% M. Vauhkonen 11.5.1994, modified from the version of J. Kaipio
% 25.4.1994. Modified 5.9.1994 by M. Vauhkonen for EIT.
% Modified 3.7.1996 by P.J. Ronkanen.
% Sparse indexing and continuous sigma integration modification by
% K. Karhunen 28.07.2006.

A = 1;


Nel=max(size(z)); %The number of electrodes.
gN=max(size(Node));
HN=max(size(Element));
M=sparse(gN,Nel);
%K=sparse(gN,gN);
s=zeros(Nel,1);
g=reshape([Node.Coordinate],3,gN)'; %Nodes
Agrad=sparse(gN,gN);
%H=reshape([Element.Topology],4,HN)';
%mH=max(max(H));

k = 1;  
Arow = zeros(4*HN,4);
Acol = zeros(4*HN,4);
Aval = zeros(4*HN,4);   

% Gauss quadrature points and weights
a=0.58541020;b=0.13819660;
ip=[b b b;a b b;b a b;b b a];

% difference matrix of the (linear) basis functions
L=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
for ii=1:HN
  % Go through all tetrahedron
  ind = Element(ii).Topology; % The indexes to g of the ii'th tetrahedron.
  gg = g(ind,:);
  ss = sigma(ind);
  int = tetraLinSigma(gg,ss,ip,L);
  id = ind(:);
  Arow(k:k+3,:) = [id id id id];
  Acol(k:k+3,:) = [ind;ind;ind;ind];
  
  if ~isempty([Element(ii).Electrode]),  % Checks if the triangle ii is the triangle that is
                                         % under the electrode.
    bind=[Element(ii).Electrode{2}]';    % Node indices under the ii'th electode    
    abc=g(bind,:);                       % Nodes under the ii'th electrode 
    InE=Element(ii).Electrode{1};        %Electrode index
    s(InE)=s(InE)+1/z(InE)*elektro([abc]);

    eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind)];
    bb1=triang2([abc]);Bb1=zeros(4,1);
    bb2=triang1([abc]);Bb2=zeros(4);

    Bb1(eind)=bb1;
    Bb2(eind,eind)=bb2;
    M(ind,InE)=M(ind,InE)-1/z(InE)*Bb1;  % minor slowdown, could be
                                         % indexed as well
    Aval(k:k+3,:) = int + 1/z(InE)*Bb2; 	
  else % The tetrahedron isn't under the electrode.
    Aval(k:k+3,:) = int; 
  end 
  k = k + 4;
end  
Agrad = sparse(Arow,Acol,Aval,gN,gN);
%Calculate next the matrix S

S=sparse(diag(s));
C=C(:,1:Nel-1);
S=C'*S*C;
M=M*C;
A=[Agrad,M;M.',S];














































