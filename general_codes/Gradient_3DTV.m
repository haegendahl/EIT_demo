% This function evaluates the gradient of the 
% approximated total variation penalty.
% July 2013
% G. Gonzalez

function grad =  Gradient_3DTV(R,theta, beta)


Nelements = size(R,1)/3;

Nnodes = length(theta);


Rx   = R(1:end/3,:);
Ry   = R(end/3+1:2*end/3,:);
Rz   = R(1+2*end/3:3*end/3,:);



grad = zeros(length(theta),1);

%elconn = cell(length(theta),1);

% for m = 1:length(theta)
%   elconn{m} = find(Rx(:,m))';
% end

for m = 1:length(theta)
    
    k = find(Rx(:,m))';
    Rxthetak = Rx(k,:)*theta;
    Rythetak = Ry(k,:)*theta;
    Rzthetak = Rz(k,:)*theta;
    grad(m) = (Rxthetak.*Rx(k,m) + Rythetak.*Ry(k,m) + Rzthetak.*Rz(k,m))'*(((Rxthetak).^2 + (Rythetak).^2 + (Rzthetak).^2 + beta).^(-.5));

end
  