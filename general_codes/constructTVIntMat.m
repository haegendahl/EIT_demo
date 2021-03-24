%%% Computes the in each element the gradient
%%% and the area of the element.
%%%
%%% Note:  Areas are computed using Heron's formula.
%%%        A stable version of the formula should be implemented
%%%        (more info in http://en.wikipedia.org/wiki/Heron's_formula) 

function [R Ai] = constructTVIntMat(g,H)

NNode=size(g,1);
NElement=size(H,1);

R = sparse(2*NElement, NNode);

Ai = zeros(NElement,1);
M = [ -1 0 1;-1 1 0];

for ii=1:NElement
    
    gg =g(H(ii,:),:);
    
      
    Adet = ((gg(3,1)-gg(1,1)) * (gg(2,2)-gg(1,2))) - (((gg(3,2)-gg(1,2)) * (gg(2,1)-gg(1,1))));
    
    Xinv = 1/Adet*[gg(2,2)-gg(1,2), gg(1,2)-gg(3,2);...
        gg(1,1)-gg(2,1) gg(3,1)-gg(1,1)];
    
      
    % Computing the area based on Heron's formula
    
   
    a = sqrt(sum((gg(2,:)-gg(1,:)).^2));
    b = sqrt(sum((gg(3,:)-gg(2,:)).^2));
    c = sqrt(sum((gg(3,:)-gg(1,:)).^2));
    
    s = (a + b + c)/2;
    
    Ai(ii) = sqrt(s*(s-a)*(s-b)*(s-c)); 


    R([ii end/2+ii],H(ii,:))=Xinv*M;
    %R(ii,H(ii,:))=sum(Xinv*M,1);
    
end

% No weighting with element size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Ai = ones(size(Ai));



