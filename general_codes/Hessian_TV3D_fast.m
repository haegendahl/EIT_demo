% This function evaluates the Hessian of the objective function
% with approximate total variation penalty.
% 
% 
 function H = Hessian_TV3D_fast(R, theta, beta)
R = sparse(R);
% 
% Rx   = R(1:end/3,:);
% Ry   = R(end/3+1:2*end/3,:);
% Rz   = R(1+2*end/3:3*end/3,:);
% 
% RxthetakALL = Rx*theta;
% RythetakALL = Ry*theta;
% RzthetakALL = Rz*theta;
% 
% coeffsALL = (RxthetakALL.^2 + RythetakALL.^2 + RzthetakALL.^2 + beta).^(-.5);
% n = numel(coeffsALL);
% D = spdiags(coeffsALL,0,n,n);
% Dx = spdiags(RxthetakALL,0,n,n);
% Dy = spdiags(RythetakALL,0,n,n);
% Dz = spdiags(RzthetakALL,0,n,n);
% RR = Dx*Rx + Dy*Ry + Dz*Rz;
% 
% fgpALL = (D*RxthetakALL)'*RxthetakALL + (D*RythetakALL)'*RythetakALL + (D*RzthetakALL)'*RzthetakALL;
% keyboard
% fpgALL = -(D.^3* RR'* RR);
% 
% H = fgpALL + fpgALL;
% 

Rx   = R(1:end/3,:);
Ry   = R(end/3+1:2*end/3,:);
Rz   = R(1+2*end/3:3*end/3,:);

Rxthetak_ALL = Rx*theta;
Rythetak_ALL = Ry*theta;
Rzthetak_ALL = Rz*theta;

tmp = (Rxthetak_ALL).^2 + (Rythetak_ALL).^2 + (Rzthetak_ALL).^2 + beta;

fgp = diag(tmp.^(-0.5))*((Rx'*Rx + Ry'*Ry + Rz'*Rz)); 
fpg = diag(-tmp.^(-1.5))*((diag(Rxthetak_ALL)*Rx)'+ (diag(Rythetak_ALL)*Ry)'+(diag(Rzthetak_ALL)*Rz)')*((diag(Rxthetak_ALL)*Rx)+ (diag(Rythetak_ALL)*Ry)+(diag(Rzthetak_ALL)*Rz));

H = fpg + fgp;

