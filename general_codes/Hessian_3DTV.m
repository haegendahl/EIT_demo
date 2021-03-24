function H = Hessian_3DTV(R, Ai, theta, beta)

Nelements = size(R,1)/3;


Rx   = R(1:end/3,:);
Ry   = R(end/3+1:2*end/3,:);
Rz   = R(1+2*end/3:3*end/3,:);


A = zeros(Nelements,1);
B = zeros(Nelements,1);
C = zeros(Nelements,1);

for k =1:Nelements
    temp1 = sqrt(((Rx(k,:)*theta).^2 + (Ry(k,:)*theta).^2+ (Rz(k,:)*theta).^2 + beta)); 
    temp2 = -sqrt(((Rx(k,:)*theta).^2 + (Ry(k,:)*theta).^2 + (Rz(k,:)*theta).^2+ beta)^3);
    
  
    temp3 = (Rx(k,:)*theta)^2;
    temp4 = (Ry(k,:)*theta)^2;
    temp5 = (Rz(k,:)*theta)^2;
    
    
    temp6 = (Rx(k,:)*theta)*(Ry(k,:)*theta);
    temp7 = (Rx(k,:)*theta)*(Rz(k,:)*theta);
    temp8 = (Ry(k,:)*theta)*(Rz(k,:)*theta);
    
    
    tt = Ai(k)/temp1; 
    
    A(k) = tt +  Ai(k)*(temp3 /temp2);
    
    B(k) =  tt +  Ai(k)*(temp4 /temp2);
    
    C(k) =  tt + Ai(k)*(temp5 / temp2);
    
    
    D(k) = (Ai(k)*temp6)/temp2;
    
    E(k) = (Ai(k)*temp7)/temp2;
    
    F(k) = (Ai(k)*temp8)/temp2;
    

    
   
    
end


H = sparse(Rx'*(diag(A)*Rx) + Ry'*(diag(B)*Ry) + Rz'*(diag(C)*Rz) + Rx'*(diag(D)*Ry) + Rx'*(diag(E)*Rz) +...
    Ry'*(diag(D)*Rx) + Ry'*(diag(F)*Rz)+ Rz'*(diag(E)*Rx) + Rz'*(diag(F)*Ry));
    




