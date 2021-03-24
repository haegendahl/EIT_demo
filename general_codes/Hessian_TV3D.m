% This function evaluates the Hessian of the objective function
% with approximate total variation penalty.


function H = Hessian_TV3D(R, theta, beta)
R = sparse(R);

Rx   = R(1:end/3,:);
Ry   = R(end/3+1:2*end/3,:);
Rz   = R(1+2*end/3:3*end/3,:);

Rxthetak_ALL = Rx*theta;
Rythetak_ALL = Ry*theta;
Rzthetak_ALL = Rz*theta;

tmp = (Rxthetak_ALL).^2 + (Rythetak_ALL).^2 + (Rzthetak_ALL).^2 + beta;

H = zeros(length(theta));
 %keyboard
for k = 1:size(Rx,1)
    Rxk=Rx(k,:);
    Ryk=Ry(k,:);
    Rzk=Rz(k,:);
    
    fgp = tmp(k)^(-.5)*((Rxk'*Rxk + Ryk'*Ryk + Rzk'*Rzk));
    
    fpg = -tmp(k)^(-1.5)*(((Rxthetak_ALL(k)*Rxk)'+ (Rythetak_ALL(k)*Ryk)'+(Rzthetak_ALL(k)*Rzk)')*((Rxthetak_ALL(k)*Rxk) + (Rythetak_ALL(k)*Ryk)+ (Rzthetak_ALL(k)*Rzk)));
       
    H = H + fpg + fgp;
end


H = sparse(H);