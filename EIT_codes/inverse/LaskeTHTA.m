function [thta,error,Gammapost,THTA] = LaskeTHTA(H,g,Element,Umeas,Gamma_n);
  
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

%%%  MESH %%%
% Next code includes bad, dynamic memory allocations. It may
% exhibit extreme slowliness on certain meshes. User discretion
% advised! KK
% :(


nH = size(H,1);
elelements = zeros(Nel,1,'uint32');
elnodes = zeros(Nel,3,'uint32');
elelements(:,1) = 2;
elnodes(:,1) = 2;

for k=1:nH
  if ~isempty(Element(k).Electrode)
    elid = Element(k).Electrode{1};
    nodeid = Element(k).Electrode{2};
    
    Ne = elelements(elid,1);    
    Nn = elnodes(elid,1);
    
    elnodes(elid,Nn:Nn+2) = nodeid;  
    elelements(elid,Ne) = k;
    elelements(elid,1) = Ne+1;    
    elnodes(elid,1) = Nn+3;
  end
end

Elind = elnodes(:,2:end);
E = elelements(:,2:end);


gN = size(g,1);
HN = size(H,1);

%%% CURRENT %%%

%[MeasPatt C] = makebhcurrentmeas;
[II1,C]=Current(16,gN,'adj'); % Measurement pattern
C = load('Tank2DOppositeInj.txt');
C=reshape(C(4:end),8,16)';

%T = zeros(Nel,8);
%cr = [1 25;3 27;5 29;7 31;9 17;11 19;13 21;15 23];
%for k=1:8, T(cr(k,:),k) = [1;-1];end

Ival = 1;
C1 = [ones(1,Nel-1);-eye(Nel-1)];
I = [zeros(gN,Ncur);C];
%I = Ival*C;

Ctde = [zeros(Nel,gN) C1];
Ihat = [I(1:gN,:);C1'*I(gN+1:gN+Nel,:)];

%%% Measurement pattern %%%

[II1,MeasPatt]=Current(16,gN,'adj'); % Measurement pattern
P = MeasPatt'*Ctde;


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

for ii = 1:Iter,
  %[A,dA_s,dA_Rez] = genmatDfast_real(g,H,Element,thta,C1);
  [A,dA_s,dA_Rez] = genmatD_real(g,H,E,Elind,thta,C1);
  PA = (A\P').';
  Uhat = PA*Ihat;
  UU = Umeas(:) - Uhat(:);
  b = A\Ihat;

  J_s = -PA*dA_s*b;  J_s = J_s(:);
  J_Rez = -PA*dA_Rez*b;  J_Rez = J_Rez(:);
  
  J = [J_s J_Rez];
  
  delta_ii = (J'*InvGamma*J)\(J'*InvGamma*UU);  
  thta = thta + (GAMMA)*delta_ii
  gammapost=inv(J'*InvGamma*J);
  Gammapost=[Gammapost gammapost];  
  THTA = [THTA thta];
  error = [error norm(UU)^2];
  figure(1),clf
  plot(Umeas(:,1)),hold on
  plot(Uhat(:,1),'r')
  drawnow
end
