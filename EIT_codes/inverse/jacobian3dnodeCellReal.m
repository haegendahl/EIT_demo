function [J] = jacobian3dnodeCellReal(arow,acol,aval,f,U,U0,style)
% Computes the Jacobian for 2D/3D EIT. Uses regrouped cell arrays
% instead of the gradient matrix Agrad (see regroupAgrad2cell).

% K. Karhunen, 09.01.2007

if nargin<6, style = 'real'; end

nn = max(size(aval));

w = 2*pi*f; 
e0 = 8.85418781762e-12; %permittivity of free space, do we really need it??

U0 = -U0.';
nm = size(U,2)*size(U0,1);
if style == 'comp'
  J = zeros(nm,2*nn);
  for ii=1:nn
    ar = arow{ii}; ac = acol{ii};
    Are = aval{ii}*U(ac,:);
    Aim = aval{ii}*U(ac+nn,:);  
    J(:,ii) = reshape(U0(:,[ar ar+nn])*[Are;Aim],nm,1);
    J(:,ii+nn) = reshape(U0(:,[ar ar+nn])*[-Aim;Are],nm,1);
  end
  J(:,1+nn:2*nn) = w*J(:,1+nn:2*nn);
elseif style == 'real'
  J = zeros(nm,nn);
  for ii=1:nn;
    JJ = U0(:,arow{ii})*(aval{ii}*U(acol{ii},:));
    JJ = JJ(:);
    J(:,ii) = JJ;
  end
end

