function [arow,acol,aval] = regroupAgrad2cell(Agrad,ng)
% Regroups the gradient matrix Agrad for the computation of
% the Jacobian in 2D/3D EIT. Regrouping tries to utilitize the full
% advantage of the matrix sparsity.

% K. Karhunen, 16.10.2006

if nargin==1
  ng = size(Agrad,2);
end

arow = cell(ng,1);
acol = cell(ng,1);
aval = cell(ng,1);
for k=1:ng
  S = reshape(Agrad(:,k),ng,ng);
  [I,J] = find(S);
  
  I = unique(I);
  J = unique(J);
  
  arow{k} = I;
  acol{k} = J;
  aval{k} = S(I,J);
end
