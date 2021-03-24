function [thta,error,Gammapost,THTA] = HomogEstimate(H,g,Element,Umeas,T,MeasPatt,Gamma_n);
% Laskee homogeenisen estimaatin johtavuudelle ja kontakti-impedanssille (reaalinen)
% Esim. Thta = HomogEstimate(H,g,Element,Uel,CurrentPattern,MeasPattern,Gamma_n)
  
maxIter = 50;

Nel = size(T,1); % elektrodien määrä
Ncur = size(T,2); % virransyöttöjen määrä

nH = size(H,1);
gN = size(g,1);
HN = size(H,1);

[II1,C1] = Current(Nel,gN,'adj'); 
C1 = C1(:,1:Nel-1); % For the voltage reference 

Ctde = sparse([zeros(Nel,gN) C1]);
II = [zeros(gN,size(T,2));C1'*T];
II1 = MeasPatt'*Ctde;


%%% INITIALS %%%

thta_0 = [.01 3]';
thta = thta_0;

GAMMA = 0.15;
THTA = [thta_0];
error = [];
Gammapost = [];

%%% RECURSION %%%

%lasketaan johtavuudet ja kontakti-impedanssit

InvGamma = inv(Gamma_n);

clear A B C D E
II1 = -II1';


ii = 1;
UU0 = Umeas;
endCriterion = 1;

while endCriterion>1e-4 & (ii<maxIter)
 
  if ii==1
    [B,C,D,E] = genmatDfast_real(g,H,Element,Nel); 
    dA_s = [B sparse(gN,Nel-1);sparse(Nel-1,gN) sparse(Nel-1,Nel-1)];
  end
  z = 1/thta(2);
  A = [thta(1)*B+z*C z*D;z*D.' z*E];
   
  if ii==1, p = symamd(A); r(p)=1:max(size(p)); end
  R = chol(A(p,p));
   
  z2 = -z^2;
  dA_z = z2*[C D;D.' E];
  
  uu = R\(R'\II1(p,:));
  MField = uu(r,:)';
  
  uu = R\(R'\II(p,:));
  Crnt = uu(r,:);
  
  Uhat = MeasPatt'*C1*uu(r(gN+1:end),:);
  UU = Umeas(:) - Uhat(:);
 
  J_s = MField*(dA_s*Crnt);
  J_z = MField*(dA_z*Crnt);
  
  J = [J_s(:) J_z(:)];  
  JIG = J'*InvGamma;
  JIGJ = JIG*J;
  
  delta_ii = (JIGJ)\(JIG*UU);  
  thta = thta + (GAMMA)*delta_ii;
  gammapost=inv(JIGJ);
  Gammapost=[Gammapost gammapost];  
  THTA = [THTA thta];
  error = [error norm(UU)^2];
  figure(1),clf
  plot(Umeas(:,1)),hold on
  plot(Uhat(:,1),'r')
  drawnow
 
  ii = ii + 1;
  endCriterion = norm(UU-UU0)/max([norm(UU) 1]);
  UU0 = UU;
  
  sz_str = sprintf(' | s = %0.2g, z = %0.2g',thta(1),thta(2));
  disp(['iteration: ' num2str(ii) ' - ' 'residual: ' num2str(endCriterion) sz_str]);
end
