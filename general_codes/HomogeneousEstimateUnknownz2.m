function [ rho0,z0,Uhat,J ] = HomogeneousEstimateUnknownz2(H,g,Node,Element,Umeas,T,MeasPatt,Gamma_n,init_sigma,init_z,REMOVE_CURRENT_INJECTING_ELECTRODE_DATA,c_inj_data_ind);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%keyboard

maxIter = 100;

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

%thta_0 = [.01 3]';
%thta_0 = [init_sigma; log(init_z)];
thta_0 = [init_sigma; init_z];
thta = thta_0;

GAMMA = 0.5;
THTA = [thta_0];
error = [];
Gammapost = [];

%%% RECURSION %%%


%keyboard
%lasketaan johtavuudet ja kontakti-impedanssit

InvGamma = inv(Gamma_n);

clear A B C D E
II1 = -II1';

ii = 1;
UU0 = Umeas;
endCriterion = 1;

while endCriterion>1e-5 & (ii<maxIter)
 
  if ii==1
    [B,C,D,E] = genmatDfast_real2D(Node,Element,C1,Nel); 
    dA_s = [B sparse(gN,Nel-1);sparse(Nel-1,gN) sparse(Nel-1,Nel-1)];
  end

  %Z = exp(thta(2));
  Z = thta(2);
  z = 1/Z;
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
  Uhat = Uhat(:);
  if REMOVE_CURRENT_INJECTING_ELECTRODE_DATA
     Uhat(c_inj_data_ind,:) = [];
  end
  UU = Umeas(:) - Uhat;
 
  J_s = MField*(dA_s*Crnt);
  J_z = MField*(dA_z*Crnt);
  J_z = J_z*exp(thta(2));
  J = [J_s(:) J_z(:)];  
  if REMOVE_CURRENT_INJECTING_ELECTRODE_DATA
     J(c_inj_data_ind,:) = [];
  end
  JIG = J'*InvGamma;
  JIGJ = JIG*J;
  
  delta_ii = (JIGJ)\(JIG*UU);  
  thta = thta + (GAMMA)*delta_ii
  %if (thta(2) < -8)
  %  thta(2) = -8;
  %end
  if (thta(2) < 1e-10)
    thta(2) = 1e-10;
  end
  
  gammapost=inv(JIGJ);
  Gammapost=[Gammapost gammapost];  
  THTA = [THTA thta];
  error = [error norm(UU)^2];
  figure(1),clf
  plot(Umeas(:,1)),hold on
  plot(Uhat(:),'r')
  drawnow
 
  ii = ii + 1;
  endCriterion = norm(UU-UU0)/max([norm(UU) 1]);
  UU0 = UU;
  
  disp(['iteration: ' num2str(ii) ' - ' 'residual: ' num2str(endCriterion)]);
end

%thta(2) = exp(thta(2));

rho0 = 1/thta(1);
z0 = thta(2);


end

