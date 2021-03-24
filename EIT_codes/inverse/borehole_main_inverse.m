% seuraava estää holtittoman käytön...
%error('Tupeloinnin esto!')

c = clock;
disp('Starting calculations')
disp(sprintf('%d.%d.%d - %d:%d:%g\n',c(3),c(2),c(1),c(4),c(5),c(6)))

%clear all

disp('Ladataan hilat...')
% $$$ load BHmesh14k Node Element H g;
% $$$ Node2=Node;Element2=Element;H2=H;g2=g;

%load boreholeInverseMesh3921.mat
load boreholeInverseMesh2057; % harva, Node, Element, H, g,
                              % priornodes priornodes2
%load boreholeInverseMesh2270 % harva hila: Node, Element, H, g
load boreholeForwardMesh14k; % tiheämpi  : Node2, Element2, H2, g2
%load boreholeForwardMesh5k

Ip = priornodes; %_coarse;
In = setdiff([1:length(g)]',Ip);


% vakioita
%gamma_kesk=0;epsilon1=0;

%load bherror003_102k2nd;  % epsilon1 gamma_kesk
%load ~/tyo06/uudet_dec2006/bherror005_randsigma.mat epsilon gamma_e

%load ~/tyo06/uudet_dec2006/bherror84k2nd003_randsigma.mat epsilon gamma_e
%load ~/tyo06/uudet_dec2006/bherror84k2nd005_randsigma.mat epsilon gamma_e
%load ~/tyo06/uudet_dec2006/bherror84k2nd005_newPrior.mat epsilon gamma_e
%load ~/tyo06/uudet_dec2006/bherror84k2nd005_newPrior5k.mat epsilon gamma_e
%load ~/tyo06/uudet_dec2006/bherror84k2nd005_newPrior_newsigma.mat epsilon gamma_e

%epsilon1 = epsilon; gamma_kesk = gamma_e;


%thta = thta(1);

%% MUUTETTU 8.12. Ok
thta = thta(1); %mean([.001 .025]);

Nel = 64;      % elektrodien määrä
Ncur = 16;     % virransyöttöjen määrä

ng1 = size(g,1);
nH1 = size(H,1);
ng2 = size(g2,1);
nH2 = size(H2,1);

maxIter = 1000;

% Johtavuuden & K-I:n odotusarvot/alkuarvaukset
sigma0 = thta;
z0 = 1*ones(Nel,1);
z = z0;

disp('Mittauskuvio & virransyöttö...')
% Mittaus- ja virtakuviot
Ival = .1; % 100 mA
[MeasPatt C] = makebhcurrentmeas;
CurrPatt = Ival*C;

%Uel = Ureference(:); 
%Uel = getReference;

%load simulaatiot_anssi/measdata_inhomogR075_dev003 Umeas;
%load simulaatiot_anssi/measurementR075dev005_smax02 Umeas;
%load simulaatiot_anssi/measurementR075dev005_smax02 Umeas;
%load forwardsim/measdata127k2nd_inhomogR075_smax02dev003.mat Umeas
%load ~/tyo06/forwardsim/measdata130k2nd_2ballsR075_dev003.mat Umeas
%load ~/tyo06/tempmeas.mat Umeas
%load ~/tyo06/uudet_dec2006/measdata102k2nd_dev005.mat Umeas
%load ~/tyo06/uudet_dec2006/measdata99k2ndR050_dev005.mat Umeas

%load ~/tyo06/uudet_dec2006/measdata92k2ndR100_dev005.mat Umeas


%load  ~/tyo06/uudet_dec2006/measdata92k2ndR100_dev005sigma10x.mat Umeas


%load ~/tyo06/uudet_dec2006/measdata96k2ndR75_eccentric_dev005.mat Umeas
 
Uel_nonoise = Umeas;




meas_noise_coef = 1e-3;      % Noise coefficient: variance of the noise is
                             %    (meas_noise_coef*(max(Uel)-min(Uel)))^2
                             %    e.g. .01 corresponds to 'one
                             %    percent noise level'

meas_noise_coef2 = 1e-2;      % Noise coefficient: variance of the noise is
                             %    var(l) = (meas_noise_coef2*Uel(l))^2
                             %    e.g. .01 corresponds to 'one percent noise level'


% $$$ [Uel_nonoise] = ComputeMeasurements(L,z,TTT,sigma_5,Node6,Element6, ...
% $$$                                     I_rms,T,MeasPatt,nM,nI,nm);


% Add gaussian noise to observations

Uel1 = Uel_nonoise + meas_noise_coef2*abs(Uel_nonoise).*randn(size(Uel_nonoise));
Uel = Uel1 + (meas_noise_coef*(max(max(Uel1))-min(min(Uel1))))*randn(size(Uel1));

%%% HUOM HUOM HUOM kerroin !!! %%%
meas_noise_coef_e = 2*meas_noise_coef; % estimated noise level (not necessarily known) before 2

beta_k = (meas_noise_coef_e*(max(max(Uel))-min(min(Uel))))^2;
var_Uel = beta_k + (meas_noise_coef2*abs(Uel)).^2;
Cv1 = var_Uel(:);
Cv = diag(Cv1(:));

disp('Kovarianssit ja Choleskyt...')
% Kohina
Gamma_n = Cv; %eye(length(Uel(:)))*1e-3;
InvGamma_n = inv(Gamma_n + gamma_kesk);
L_n = chol(InvGamma_n);

figure(1), clf, plot(sqrt(diag(Gamma_n))), drawnow
figure(2), clf, plot(sqrt(diag(gamma_kesk))), drawnow
figure(3), clf, plot(epsilon1), drawnow

% priori
%R = 5*MakeGradientRegMat3DFast(Node,Element);
%InvGamma_pr = R'*R;
%L_pr = R;%chol(InvGamma_pr);

alpha = 800;
cond = sigma0;
range = [.001 .025];
[eta0 Igamma] = ConstructProperSmoothnessPrior3D(g,H,Node,Element,cond,alpha,range,Ip);
s = diag(inv(Igamma)); [foo I] = sort(s(In)); In2 = In(I);
%figure(4), plot(s([Ip;In2]))
figure(3), plot(In,s(In),'b.',Ip,s(Ip),'r.')
drawnow
%return

L_pr = chol(Igamma);

%[InvGamma_pr Gamma_pr] = ConstructBHSmoothnessPrior3D(Element,Node,g,H,thta_real,priornodes); 
%L_pr=chol(InvGamma_pr);

%[InvGamma_pr mu_pr Gamma_pr] = ConstructSmoothnessPrior3Dgeo(Element,Node,g,H,thta_real); 
%L_pr=chol(InvGamma_pr);
% w = chol(Gamma_pr)'*randn(N,1);  random draw from the smoothness prior

%%%%%%%%%%%%%%%%%%%%%
%% INVERSE PROBLEM %%
%%%%%%%%%%%%%%%%%%%%%

disp('Interpolaatiomatriisi...')
%sigma1 = eta0; %
sigma1 = sigma0*ones(ng1,1);

%sigma2 = sigma0*ones(ng2,1);

WW = speye(ng1,ng2);

%load meshmap2k_to_14k P2

[ff P2 Emap] = Interpolate2Newmesh3D(g,H,sigma1,g2,[],[]);
%load meshmap2k_to_14k P2

sigma2 = P2*sigma1;

disp('Suora ongelma...')
% Lineaarinen kanta

%load homogvoltagedata0029 Uref2
Uref2 = ForwardSolution3d(Node2,Element2,CurrPatt,[],sigma2,z,MeasPatt,'real');
%load homogvoltagedata0029_bh5k Uref2 

Current = Uref2.Current;
MsField = Uref2.MeasField;

disp('Jacobin laskenta ja projisointi...')
% Jacobi lineaarisessa kannassa

WW = speye(ng2);
Agrad = jacobian3dnodeFast(Node2,Element2,WW);
%load JacobianAgrad_eta12102006_mesh2057.mat Agrad

%load jacobiandata2057sigma0029.mat ar ac av
[ar ac av] = regroupAgrad2cell(Agrad,ng2);
%load jacobiandata2057sigma0029_bh5k.mat ar ac av

J2 = jacobian3dnodeCellReal(ar,ac,av,[],Current,MsField,'real');
%J2 = jacobian3dnodeFast(Node2,Element2,WW,Agrad,Current,MsField,'real');

%load JacobianAgrad_eta12102006_mesh2057.mat Agrad J2
%load JacobianAgrad_eta09102006_mesh1474.mat Agrad J2
%load JacobianAgrad_eta0_mesh1474.mat Agrad J2
%load JacobianAgrad_sigma0021_mesh1474.mat Agrad J2
%load JacobianAgrad_sigma0021_mesh3921.mat Agrad J2
%load JacobianAgrad_sigma0021_mesh2270.mat Agrad J2
%load JacobianAgradR075_mesh1804.mat Agrad J2
disp('Jacobian (dense) and Agrad loaded!')
%save JacobianAgrad_sigma0021_mesh2270.mat Agrad J2
%save JacobianAgrad_sigma0021_mesh3921.mat Agrad J2
%disp('Jacobian (dense) and Agrad saved!')

J = J2*P2;

sigma = zeros(ng1,1);
sigma(:,1) = sigma1;

Urefel = Uref2.Electrode(:);

Fnorm(1) = norm(L_n*(Uel-Urefel-epsilon1))^2 + norm(L_pr*(sigma(:,1)-sigma1))^2;
measnorm(1) = norm(Uel-Urefel-epsilon1);

args = {Node2,Element2,CurrPatt,z,MeasPatt,L_n,L_pr,P2};

disp('Käänteisongelman iteratiivinen ratkaisu...')
a = .5;

lsearchsteps = 8;
%sigma_err = zeros(maxIter,1);

lminlast = 1;
for kk=2:maxIter
  
   HH=[L_n*J;L_pr];
   zz=[L_n*(Uel-Urefel-epsilon1);-L_pr*(sigma(:,kk-1)-sigma1)];
   suunta = HH\zz;
   
%   sigma(:,kk) = sigma(:,kk-1) + a*suunta;
%   step = 2;
   
   disp('linesearch...');
   [beta jj] = sigmastep(sigma(:,kk-1),suunta,.0005);
   
   bmax = max(beta);

   fsteplength = bmax/lsearchsteps;    
   fstepvec = zeros(1,lsearchsteps);
   fstepvec(1:3) = linspace(min([fsteplength/3 1e-5]),fsteplength,3);
   fstepvec(3:end) = linspace(fstepvec(3),bmax,lsearchsteps-2);
   
   lsearch_go = 1;
   fstep = 0; fsall = 0;
   cnt = 0;  
   F = 0;
   while lsearch_go & (cnt<lsearchsteps)
     cnt = cnt + 1;
     disp(cnt)
     
     fstep = fstepvec(cnt);
     fsall(cnt) = fstep;
     ss1 = sigma(:,kk-1) + (beta*fstep/bmax).*suunta;
     F(cnt) = LaskeFnormi(ss1,Uel-epsilon1,sigma1,sigma(:,kk-1),args);     
     [ff fI] = sort(F(1:cnt));
     if (cnt>2) & ((fI(1)~=1) & (fI(1)~=cnt)), lsearch_go = 0; end     
     %if (fstep+fsteplength*1.5)>bmax, fsteplength = bmax-fstep; end     
   end
   if (fI(1)==1) | (fI(1)==cnt) 
     lmin = fstepvec(fI(1)); 
   else, 
     p = fI(1);
     p = polyfit(fstepvec(p-1:p+1),F(p-1:p+1),2);
     lmin = -.5*p(2)/p(1);           
   end   

%%%%%%%%%%%
%lmin = 0.5;
%%%%%%%%%%%

   bstep(:,kk) = beta*lmin/bmax;
  
   
   disp(['ASKEL: ' num2str(max(bstep(:,kk)))])
   disp(['ZERO : ' num2str(length(find(beta==0)))])
   
   sigma(:,kk) = sigma(:,kk-1) + bstep(:,kk).*suunta;
   
% $$$    I = find(sigma(:,kk) <= 0);
% $$$    if ~isempty(I), disp(['NEGATIVE or ZERO: ' num2str(length(I)) ...
% $$$                          ' max sigma: ' num2str(max(sigma(:,kk)))]), end
% $$$    sigma(I,kk) = 0.001;
   sigma2 = P2*sigma(:,kk);
 
   sigma_err(kk) = norm(sigma(:,kk)-sigma(:,kk-1))/max([1 norm(sigma(:,kk))]);
   suunta_err(kk) = norm(bstep(:,kk).*suunta);
   disp([sigma_err(kk) suunta_err(kk)])
   
   figure(1),clf
   plot(fsall,F,'b-+'), hold on
   v = axis;
   plot([lmin lmin],v(3:4),'r-')
   drawnow
   
   figure(2), clf
   plot(sigma(:,kk))
   drawnow
   
   figure(3), clf
   BH_sliceplot(g,sigma(:,kk));
   drawnow
   
   save ~/tyo06/sigmaIter.mat sigma
   
   %if (bmax>1e-5) & (lmin<=1.1e-5)
   if(suunta_err(kk) < 1e-6)
     disp('Break!')
     break
   end
   
   if kk<maxIter
     Uref2 = ForwardSolution3d(Node2,Element2,CurrPatt,[],sigma2,z,MeasPatt,'real');
     J2 = jacobian3dnodeCellReal(ar,ac,av,[],Uref2.Current,Uref2.MeasField,'real');
%     J2 = jacobian3dnodeFast(Node2,Element2,WW,Agrad,Uref2.Current, ...
%                             Uref2.MeasField,'real');
     
     Urefel = Uref2.Electrode(:);

     Fnorm(kk) = norm(L_n*(Uel-Urefel-epsilon1))^2 + norm(L_pr*(sigma(:,kk)-sigma1))^2;
     measnorm(kk) = norm(Uel-Urefel-epsilon1);
     figure(1), clf
     subplot(2,1,1),plot(fsall,F,'b-+'), v=axis; hold on,plot([lmin lmin],v(3:4),'r-'),subplot(2,1,2),plot(Fnorm)
     figure(2), clf
     subplot(2,1,1), plot(measnorm), subplot(2,1,2), plot(1:640,Uel,'b-',1:640,Urefel+epsilon1,'r-')
     drawnow
     
     J = J2*P2;
   end
   
   lminlast = lmin;
   disp(sprintf('Iteration: %d/%d done!\n',kk,maxIter))
end
sigma_err = sigma_err(2:kk);
suunta_err = suunta_err(2:kk);

disp('Saving...')
%save ~/tyo06/testiInverse.mat sigma -v6
save ~/tyo06/sigma_stand.mat sigma sigma_err -v6

c = clock;
disp('All done!')
disp(sprintf('%d.%d.%d - %d:%d:%g\n',c(3),c(2),c(1),c(4),c(5),c(6)))

return
