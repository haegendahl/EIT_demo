function [J] = jacobian3dnodeRealCell(arow,acol,aval,WW,U,U0);
% Computes the Jacobian for 3D EIT. Uses regrouped cell arrays
% instead of the gradient matrix Agrad.

% K. Karhunen, 16.10.2006

nn = max(size(aval));
J = zeros(size(U,2)*size(U0,2),nn);

U0 = -U0.';
for ii=1:nn;
  JJ = U0(:,arow{ii})*(aval{ii}*U(acol{ii},:));
  JJ = JJ(:);
  J(:,ii) = JJ;
end
