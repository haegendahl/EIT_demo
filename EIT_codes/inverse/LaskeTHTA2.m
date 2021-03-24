function [thta,error,Gammapost,THTA] = LaskeTHTA2(H,g,Element,Umeas,Gamma_n);
  
  % Ohjelma laskee kohteen ominaisjohtokyvyn ja suhteellisen
  % permittiivisyyden seka kontakti-impedanssin Gauss-Newton
  % menetelmalla.  
  %
  % T. Vilhunen 5.3.2001 
  % Modified by LH 12.08.2003  
  % LH 26.11.2003
  % Last modification, 20.1.04 LH
  
thta = 1;
Nel = 16; % elektrodien määrä
Ncur = 8; % virransyöttöjen määrä

nH = size(H,1);
gN = size(g,1);
HN = size(H,1);

%%% CURRENT %%%

T = load('Tank2DOppositeInj.txt');
T = reshape(T(4:end),8,16)';  % currentpattern

[II1,C1] = Current(16,gN,'adj'); 
C1 = C1(:,1:Nel-1); % For the voltage reference 

Ctde = sparse([zeros(Nel,gN) C1]);
II = [zeros(gN,size(T,2));C1'*T];


%%% Measurement pattern %%%
MeasPatt = toeplitz([1;-1;zeros(14,1)],[1 zeros(1,14) -1]);
II1 = MeasPatt'*Ctde;


%%% INITIALS %%%

thta_0 = [1 1]';
thta = thta_0;

GAMMA = 0.5;
THTA = [thta_0];
error = [];
Gammapost = [];

%%% RECURSION %%%

%lasketaan johtavuudet ja kontakti-impedanssit

Iter = 20;

InvGamma = inv(Gamma_n);

clear A B C D E
II1 = -II1';

for ii = 1:Iter,
% [A,dA_s,dA_Rez] = genmatD_real(g,H,E,Elind,thta,C1);
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
  thta = thta + (GAMMA)*delta_ii
  gammapost=inv(JIGJ);
  Gammapost=[Gammapost gammapost];  
  THTA = [THTA thta];
  error = [error norm(UU)^2];
  figure(1),clf
  plot(Umeas(:,1)),hold on
  plot(Uhat(:,1),'r')
  drawnow
end
