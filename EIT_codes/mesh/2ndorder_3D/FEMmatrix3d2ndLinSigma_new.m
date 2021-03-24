function [AA]=FEMmatrix3d2ndLinSigma(g,H,Eind,elind,z,sigma);

%keyboard

% FOR SECOND ORDER BASIS FUNCTIONS
% function [AA]=genmat2quad(g,gp,Hp,E,rho,z) supplies matrices
% needed for the FEM solution in the EIT forward problem. 
% g contains the nodes in (x,y) coordinates, gp in polar coordinates, 
% Hp the indexes to g of the triangle vertices.
% E contains the indices of the triangles that have
% three nodes on the boundary and are under the electrodes. 
% Needed for calculations of boundary integrals. rho is
% the resistivity vector. Einds tells which nodes are on the 
% boundary.
% M. Vauhkonen 11.5.1994, modified from the version of J. Kaipio
% 25.4.1994. Modified 5.9.1994 by M. Vauhkonen for EIT.
% Modified 31.5.2000 by P.J. Vauhkonen for 3D EIT. 
aa=fieldnames(Eind);
bb=fieldnames(elind);
Nel=max(size(aa)); %The number of electrodes.
gN=max(size(g));
HN=max(size(H));
M=sparse(gN,Nel);
K=sparse(gN,gN);
S=sparse(Nel,Nel);
no=10;
ne=6; 
nn=size(H,2);
A=sparse(gN,gN); 

for ii=1:HN
  %A=sparse(gN,gN); 
  % Go through all triangles
  ind=H(ii,:); % The indices to g of the ii'th triangle.
  gg=g(ind,:);   % The 3x2 matrix of triangle vertices in (x,y) coord.
  tedratot=tedraLinSigma(gg,sigma(ind(1:4)));
  fE1=[];
  fE2=[];
  for kk=1:Nel
   pp=find(eval(['Eind.' aa{kk}])==ii); %Checks if the tetrahedron ii 
    if ~isempty(pp)                     %is the triangle that is 
     fE1=kk;fE2=pp;                             % under the electrode. 
    end
  end
  if ~isempty(fE2),
  el=reshape(eval(['elind.' bb{fE1}]),ne,size(eval(['elind.' bb{fE1}]),2)/ne)';
  gel=el(fE2,:);
  %keyboard
  %S(fE1,fE1)=S(fE1,fE1)+1/z(fE1)*ele(g(gel,:));
  S(fE1,fE1)=S(fE1,fE1)+1/z(fE1)*elektro(g(gel(1:3),:));
  %elektro(g(gel(1:3),:))
  %ele(g(gel,:))
  bb1=[];
  for i=1:length(gel)
     bb1=[bb1;find(ind==gel(i))];
  end
  Bb1=zeros(nn,1);
  Bb2=zeros(nn);
  Bb1(bb1)=tri21(g(gel,:));
  Bb2(bb1,bb1)=tri22(g(gel,:));
    %keyboard
    M(ind,fE1)=M(ind,fE1)-1/z(fE1)*Bb1;
    K(ind,ind)=K(ind,ind)+1/z(fE1)*Bb2;
    A(ind,ind)=A(ind,ind)+tedratot;
  else %The triangle isn't under the electrode.
    A(ind,ind)=A(ind,ind)+tedratot;
  end
end  
%Calculate next the matrix 

%%
[I,C]=current_eit(Nel,1,'cee',g);
clear I
C=C(:,1:Nel-1);
S=C'*S*C;
M=M*C;
AA=sparse([K+A,M;M',S]);

% Symmetric completions

%K=K+tril(K,-1)';
%A=[K,M;M',S];
