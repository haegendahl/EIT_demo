% This function evaluates the Hessian of the 
% approximated total variation penalty.
% July 2013
% G. Gonzalez

 function H = Hessian_TV3D_fast_2(R, theta, beta)
R = sparse(R);

Rx   = R(1:end/3,:);
Ry   = R(end/3+1:2*end/3,:);
Rz   = R(1+2*end/3:3*end/3,:);

Rxtht = Rx*theta;
Rytht = Ry*theta;
Rztht = Rz*theta;

tmp = ((Rxtht).^2 + (Rytht).^2 + (Rztht).^2 + beta);

n = numel(tmp);

tmp_sp1 = spdiags(tmp.^(-.5),0,n,n);

a = spdiags(-tmp.^(-1.5).*Rxtht,0,n,n);
b = spdiags(-tmp.^(-1.5).*Rytht,0,n,n);
c = spdiags(-tmp.^(-1.5).*Rztht,0,n,n);

t = spdiags(Rxtht,0,n,n);
u = spdiags(Rytht,0,n,n);
v = spdiags(Rztht,0,n,n);

fgp = (tmp_sp1*Rx)'*Rx + (tmp_sp1*Ry)'*Ry + (tmp_sp1*Rz)'*Rz; 

fpg = ((a*Rx)'+ (b*Ry)'+(c*Rz)')*(t*Rx + u*Ry + v*Rz);

H = fpg + fgp;


