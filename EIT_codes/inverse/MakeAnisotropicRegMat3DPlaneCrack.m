function R = MakeAnisotropicRegMat3DPlaneCrack(Node,Element,glim,lambda)

R = 1;

% Function R=regmatnode(Node,Element);
% calculates a regularization matrix in the
% case where the conductivity distribution is
% estimated in the nodal basis.

% Aku Sepp?nen 30.8.2002

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
   Qflag = [0 0 0];
   
   if (sum(abs(gg(:,1)-glim(1,1))<1e-4)==3) | (sum(abs(gg(:,1)-glim(1,2))<1e-4)==3)
     Q(2,:) = lambda(2)*Q(2,:);  % y-direction
     Q(3,:) = lambda(3)*Q(3,:);  % z-direction
     Qflag([2 3]) = 1;
   end
   
   if (sum(abs(gg(:,2)-glim(2,1))<1e-4)==3) | (sum(abs(gg(:,2)-glim(2,2))<1e-4)==3)
     Q(1,:) = lambda(1)*Q(1,:);  % x-direction
     if ~Qflag(3), Q(3,:) = lambda(3)*Q(3,:); end % z-direction
     Qflag([1 3]) = 1;
   end
   
   cp = mean(gg(:,1:2));
   if (sum(abs(gg(:,3)-glim(3,2))<1e-4)==3)
     if ~((abs(25-cp(1))<5) & (abs(25-cp(2))<11))
       if ~Qflag(1), Q(1,:) = lambda(1)*Q(1,:); end % x-direction
       if ~Qflag(2), Q(2,:) = lambda(2)*Q(2,:); end % y-direction
       Qflag([1 2]) = 1;       
     else   
       Q(1,:) = 1e-2*Q(1,:);
       Q(2,:) = Q(2,:);       
     end
     
   end
   
   if (sum(abs(gg(:,3)-glim(3,1))<1e-4)==3)
       if ~Qflag(1), Q(1,:) = 0.1*lambda(1)*Q(1,:); end % x-direction
       if ~Qflag(2), Q(2,:) = 0.1*lambda(2)*Q(2,:); end % y-direction
       Qflag([1 2]) = 1;
   end
   
   if ~any(Qflag)
     Q(3,:) = 1*Q(3,:);  %z-direction
   end
   
   if ((abs(25-cp(1))<5) & (abs(25-cp(2))<11))
     Q(1,:) = 1e-2*Q(1,:); 
     Q(2,:) = 1*Q(2,:); 
     Q(3,:) = 1e-2*Q(3,:);
   end
      
   
   Rcol(k:k+2,:) = [elind;elind;elind];
   Rval(k:k+2,:) = Q;       % ei painotusta pinta-alalla!!!!!!!!!
   k=k+3;
end
Rrow = repmat([1:3*NElement]',1,4);
R=sparse(Rrow,Rcol,Rval,3*NElement,NNode);
