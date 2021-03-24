function [A,dA_s,dA_Rez] = genmatD_real(gn,Hn,Eind,elind,thta,C); 
  
A =1; 
dA_s = 1; 
dA_Rez = 1;

  % Ohjelma laskee fem-matriisin ja sen derivaatat estimoitavien
  % parametrien suhteen.
  %
  % T. Vilhunen 5.3.2001

%j=sqrt(-1);

%eps_r = thta(1);
sigma = thta(1);
Rez = thta(2);
%Imz = thta(4);

Nel = size(elind,1);
gN = max(size(gn));
HN = max(size(Hn));

K = sparse(gN,gN);
M = sparse(gN,Nel);
s = 0;

%K_e = sparse(gN,gN);
%M_e = sparse(gN,Nel);
K_s = sparse(gN,gN);
M_s = sparse(gN,Nel);
K_Rez = sparse(gN,gN);
M_Rez = sparse(gN,Nel);
%K_Imz = sparse(gN,gN);
%M_Imz = sparse(gN,Nel);
%s_e = 0;
s_s = 0;
s_Rez = 0;
%s_Imz = 0;
acnt=0;
for ii = 1:HN,
   ind = (Hn(ii,:));
   gg = gn(ind,:);
   %tedratot = (sigma + j*w*eps_0*eps_r)*tedra(gg);
   tedratot = (sigma)*tedra(gg);
   %tedratot_e = j*w*eps_0*tedra(gg);
   tedratot_s = tedra(gg);
   [fE1,fE2] = find(Eind==ii);
   fE = [fE1,fE2];
   if ~isempty(fE),
      el = reshape(elind(fE1,:),3,size(elind,2)/3)';
      gel = el(fE2,:);
      if fE1==1
         %s = s + 1/(Rez+j*Imz)*elektro(gn(gel,:));
         s = s + 1/(Rez)*elektro(gn(gel,:));
         %s_e = s_e + 0;
         s_s = s_s + 0;
         s_Rez = s_Rez - 1/(Rez)^2*elektro(gn(gel,:));
         %s_Rez = s_Rez - 1/(Rez+j*Imz)^2*elektro(gn(gel,:));
         %s_Imz = s_Imz - j/(Rez+j*Imz)^2*elektro(gn(gel,:));
      end
      for il = 1:4,
         IL = find(ind(il)==gel);
         if ~isempty(IL),
            %M(ind(il),fE1) = M(ind(il),fE1) - 1/(Rez+j*Imz)* ...
            %   triang22_old(IL,gn(gel,:));
            M(ind(il),fE1) = M(ind(il),fE1) - 1/(Rez)*triang22_old(IL,gn(gel,:));            
            %M_e(ind(il),fE1) = M_e(ind(il),fE1) + 0;
            M_s(ind(il),fE1) = M_s(ind(il),fE1) + 0;
            %M_Rez(ind(il),fE1) = M_Rez(ind(il),fE1) + 1/(Rez+j*Imz)^2* ...
            %   triang22_old(IL,gn(gel,:));
            M_Rez(ind(il),fE1) = M_Rez(ind(il),fE1) + 1/(Rez)^2* ...
               triang22_old(IL,gn(gel,:));
            %M_Imz(ind(il),fE1) = M_Imz(ind(il),fE1) + j/(Rez+j*Imz)^2* ...
            %   triang22_old(IL,gn(gel,:));
         end
         for im = 1:4,
            JL = find(ind(im)==gel);
            if ~isempty(IL) & ~isempty(JL),
               %K(ind(il),ind(im)) = K(ind(il),ind(im)) + tedratot(il,im) + ...
               %   1/(Rez+j*Imz)*triang1_old(IL,JL,gn(gel,:));
               K(ind(il),ind(im)) = K(ind(il),ind(im)) + tedratot(il,im) + ...
                  1/(Rez)*triang1_old(IL,JL,gn(gel,:));
               %K_e(ind(il),ind(im)) = K_e(ind(il),ind(im)) + tedratot_e(il,im);
               K_s(ind(il),ind(im)) = K_s(ind(il),ind(im)) + tedratot_s(il,im);
               %K_Rez(ind(il),ind(im)) = K_Rez(ind(il),ind(im)) - ...
               %   1/(Rez+j*Imz)^2*triang1_old(IL,JL,gn(gel,:));
               K_Rez(ind(il),ind(im)) = K_Rez(ind(il),ind(im)) - ...
                  1/(Rez)^2*triang1_old(IL,JL,gn(gel,:));
               %K_Imz(ind(il),ind(im)) = K_Imz(ind(il),ind(im)) - ...
               %   j/(Rez+j*Imz)^2*triang1_old(IL,JL,gn(gel,:));
            else
               acnt = acnt + 1;
               K(ind(il),ind(im)) = K(ind(il),ind(im)) + tedratot(il,im);
               %K_e(ind(il),ind(im)) = K_e(ind(il),ind(im)) + tedratot_e(il,im);
               K_s(ind(il),ind(im)) = K_s(ind(il),ind(im)) + tedratot_s(il,im);
               K_Rez(ind(il),ind(im)) = K_Rez(ind(il),ind(im)) + 0;
               %K_Imz(ind(il),ind(im)) = K_Imz(ind(il),ind(im))+0;
            end
         end
      end
   else 
      for il = 1:4,
         for im = 1:4,
            K(ind(il),ind(im)) = K(ind(il),ind(im)) + tedratot(il,im);
            %K_e(ind(il),ind(im)) = K_e(ind(il),ind(im)) + tedratot_e(il,im);
            K_s(ind(il),ind(im)) = K_s(ind(il),ind(im)) + tedratot_s(il,im);
            K_Rez(ind(il),ind(im)) = K_Rez(ind(il),ind(im)) + 0;
            %K_Imz(ind(il),ind(im)) = K_Imz(ind(il),ind(im)) + 0;
         end
      end
   end
end
disp(acnt)

S = sparse(diag(s*ones(Nel,1)));
%S_e = sparse(diag(s_e*ones(Nel,1)));
S_s = sparse(diag(s_s*ones(Nel,1)));
S_Rez = sparse(diag(s_Rez*ones(Nel,1)));
%S_Imz = sparse(diag(s_Imz*ones(Nel,1)));

S = C'*S*C;
%S_e = C'*S_e*C;
S_s = C'*S_s*C;
S_Rez = C'*S_Rez*C;
%S_Imz = C'*S_Imz*C;
M = M*C;
%M_e = M_e*C;
M_s = M_s*C;
M_Rez = M_Rez*C;
%M_Imz = M_Imz*C;

A = sparse([K,M;M.',S]);
%dA_e = sparse([K_e,M_e;M_e.',S_e]);
dA_s = sparse([K_s,M_s;M_s.',S_s]);
dA_Rez = sparse([K_Rez,M_Rez;M_Rez.',S_Rez]);
%dA_Imz = sparse([K_Imz,M_Imz;M_Imz.',S_Imz]);
