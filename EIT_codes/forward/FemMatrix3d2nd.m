function [A]=FemMatrix3d2nd(Node,Element,sigma,z,C,mtype);

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
nn=10; 

%H=reshape([Element.Topology],4,HN)';
%mH=max(max(H));

k = 1;
Arow = zeros(10*HN,10);
Acol = zeros(10*HN,10);
Aval = zeros(10*HN,10);
%% 4 point integration (2nd order accurate) %%
wght=[1/24*ones(4,1)];
a=0.58541020;b=0.13819660;
%intp=[b b b;a b b;b a b;b b a];
intp = [a b b;b a b;b b a;b b b];
%%


if strcmp(upper(mtype),'NETGEN')
  %disp(['** Using NETGEN-type ordering of element nodes (see sfl-vector in
  %FemMatrix3d2ndLinSigma.m)']) 
  sfl = [1:5 8 6 7 9 10];  % NETGEN only, the colums of H matrix must be applied in this order for the
                         % integration
else
  sfl = [1:10];
end


np = size(intp,1);
for ii=1:HN
  % Go through all tetrahedron
  ind = Element(ii).Topology; % The indexes to g of the ii'th tetrahedron.
  ind = ind(sfl);
  gg = g(ind,:);
  ss = sigma(ind(:));
  int = tetra2nd(gg,ss,wght,intp,np); 
  Acol(k:k+9,:) = [ind;ind;ind;ind;ind;ind;ind;ind;ind;ind];ind=ind.';
  Arow(k:k+9,:) = [ind ind ind ind ind ind ind ind ind ind];  
  
  if ~isempty([Element(ii).Electrode]),  % Checks if the triangle ii is the triangle that is
                                         % under the electrode.
    bind=[Element(ii).Electrode{2}];    % Node indices under the
                                        % ii'th electode    
%    bind = bind([1 4 2 5 3 6]);
    abc=g(bind,:);       % Nodes under the ii'th electrode 
    InE=Element(ii).Electrode{1};        %Electrode index
    s(InE)=s(InE)+1/z(InE)*elektro2nd(abc);

    eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind), ...
          find(bind(4)==ind),find(bind(5)==ind),find(bind(6)==ind)]; 
    bb1=triang2nd1([abc]);Bb1=zeros(nn,1);
    bb2=triang2nd2([abc]);Bb2=zeros(nn);
   
    Bb1(eind)=bb1;
    Bb2(eind,eind)=bb2;
    M(ind,InE)=M(ind,InE)-1/z(InE)*Bb1;  % minor slowdown, could be
                                         % indexed as well
    Aval(k:k+9,:) = int + 1/z(InE)*Bb2; 	
  else % The tetrahedron isn't under the electrode.
    Aval(k:k+9,:) = int; 
  end 
  k = k + 10;
end  
Agrad = sparse(Arow,Acol,Aval,gN,gN);
%Calculate next the matrix S

S=sparse(diag(s));
C=C(:,1:Nel-1);
S=C'*S*C;
M=M*C;
A=[Agrad,M;M.',S];

