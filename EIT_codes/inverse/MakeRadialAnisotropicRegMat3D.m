function R = MakeRadialAnisotropicRegMat3D(Node,Element,lambda,radius)

R = 1;

% Function R=regmatnode(Node,Element);
% calculates a regularization matrix in the
% case where the conductivity distribution is
% estimated in the nodal basis.

% Aku Seppänen 30.8.2002

% Sparse indexing update to speed up the calculation, K. Karhunen
% 27.07.2006 

NNode=max(size(Node));
NElement=max(size(Element));
g=reshape([Node.Coordinate],3,NNode)';
M=[-1 0 0 1;-1 0 1 0;-1 1 0 0];
k=1;

Rcol = zeros(3*NElement,4);
Rval = zeros(3*NElement,4);
Lambda = diag(lambda);

for ii=1:NElement
   elind = Element(ii).Topology; % Indices of the nodes.
   gg = g(elind,:);
   
   %V = (1/6)*det([ones(4,1) gg]);
   X=[gg(4,1)-gg(1,1) gg(4,2)-gg(1,2) gg(4,3)-gg(1,3); ...
      gg(3,1)-gg(1,1) gg(3,2)-gg(1,2) gg(3,3)-gg(1,3); ...
      gg(2,1)-gg(1,1) gg(2,2)-gg(1,2) gg(2,3)-gg(1,3)];
   
   Q=inv(X)*M;
   %R(k:k+2,elind)=V*Q;
   
   r = sqrt(sum(gg(:,1:2).^2,2));
   cp = mean(gg);  % center of the tetrahedron
  [TH RR] = cart2pol(cp(1),cp(2));
   if any(abs(r-radius)<1e-2)
     if (TH<0) | (TH>pi/6)
       v1 = cp(:)/norm(cp);
       v2 = [-v1(2);v1(1);0];
       v3 = [0;0;1];
       V = [v1 v2 v3];
       W = V*Lambda*V';
       Q = W*Q;
     end     
   else     
     Q(1,:) = .5*Q(1,:);  % x-direction
     Q(2,:) = .5*Q(2,:);  % y-direction
     Q(3,:) = .5*Q(3,:);  % z-direction
   end

   if (TH>0) && (TH<pi/6)
       v1 = cp(:)/norm(cp);
       v2 = [-v1(2);v1(1);0];
       v3 = [0;0;1];
       V = [v1 v2 v3];
       W = V*diag([100 .1 50])*V';
       Q = W*Q;    
   end
   
   Rcol(k:k+2,:) = [elind;elind;elind];
   Rval(k:k+2,:) = Q;       % ei painotusta pinta-alalla!!!!!!!!!
   k=k+3;
end
Rrow = repmat([1:3*NElement]',1,4);
R=sparse(Rrow,Rcol,Rval,3*NElement,NNode);
