
function [sigma_z,Uhat]=Estimate_homo_contact(g1st,H1st,Node1st,Element1st,zmax_coef,z_bg_max_coef,Nel,CurrentPattern,MeasPattern,eltetra1st,E1st,InvGamma_n,Uel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Estimate contact impedance          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Estimate the initial guess for the whole est of conductivity
% Treat as homo. case


sN=max(size(g1st));
Ncur=32;

%%%%  Current pattern %%%%

Ival = 1; % virrat tarkast ett� samat muissakin filess�
T=CurrentPattern;



%%% Load covairance information of noise



% Uel=mit11(:,20);
% gamma_e=diag(gamma_e);
% gamma_e=diag(gamma_e);
%
% InvGamma_n=inv(gamma_e);



%%%Expected value
sigma_0=1.5;
z_0=0.1;


WW=speye(sN,sN);
Agrad=jacobian3dnode(Node1st,Element1st,WW);  %By Antti

%%%% Agrad=jacobian3dnodeFast(Node,Element,WW); % Change by Dong

ng = max(size(Node1st));
[ar,ac,av] = regroupAgrad2cell(Agrad,ng);

iter=30;

% initial guess
z_intial=0.5
z=z_intial*ones(Nel,1);
sigma=1.3;

k=0.5; % askel
step1=0.1;
bet=1000; % t�m� simuloiduilla mittauksilla 1000  Tokhonove reg para...
kriteeri=1;


% -- Contact impedance prior

z0 = 0;
zmax = zmax_coef*z_intial;
zmin = z0 - zmax;
zrange = [zmin zmax];
Gamma_z = ((diff(zrange)/6)^2) * eye(Nel);
% background

zbg0 = z_intial;
zbgmax = z_bg_max_coef*z_intial;
zbgmin = 2*zbg0 - zbgmax;
zbgrange = [zbgmin zbgmax];
bg_var = (diff(zbgrange)/6)^2;
Gamma_z = Gamma_z + bg_var*ones(Nel);
IGamma_z = inv(Gamma_z);

%%%%

thta=[sigma z'];
thta=thta(:);
beta=log(thta);


z_0=z_0*ones(Nel,1);
thta_0=[sigma_0 z_0'];
thta_0=thta_0(:);
beta_0=log(thta_0);



z_ii=zeros(Nel,iter);
sigma_ii=zeros(iter,1);
virhe=zeros(iter,1);

tic

%for ii=1:iter,
ii=1;

IGamma_sigma_z=eye(Nel+1);

IGamma_sigma_z(2:Nel+1,2:Nel+1)=IGamma_z ;

while (kriteeri>(10^(-16)) & (ii<iter))
    ii
    sss=sigma*ones(max(size(Node1st)),1);
    %     Uref=ForwardSolution3dnode_KIAN(Node,Element,T,sss,z,MeasPatt,'real');
    %     %% my�s nopeutettu versio  antti
    %    Uref=ForwardSolution3dKIAN_2(Node1st,Element1st,T,sigma,sss,z,MeasPattern,'real'); %% my�s nopeutettu versio Dong
    Uref = ForwardSolution3dCI(Node1st,Element1st,CurrentPattern,sss,z,MeasPattern,'real');
    %     virhe(ii)=norm(Uel-Uref.Electrode(:))^2;
    zz=cell(3,1);
    zz{1}=Uref.Electrode;
    zz{2}=Uref.Current;
    zz{3}=Uref.MeasField;
    %     J_z=ContactImpedanceJacobi(g1st,H1st,z,zz,E1st,eltetra1st);
    J_z = ContactImpedanceJacobi(g1st,H1st,z,zz,eltetra1st,E1st);
    J=jacobian3dnodeRealCell(ar,ac,av,WW,Uref.Current(1:sN,:),Uref.MeasField(1:sN,:));
    Jtilde=[sum(J,2),J_z];
    %     Gamma_post=inv(Jtilde'*InvGamma_n*Jtilde+bet*eye(Nel+1));
    
    % logaritminen
    
    Jtilde=Jtilde*diag(exp(beta));
    
    
    % tikhonov regularisoitu
    
    %deltatheta=(Jtilde'*InvGamma_n*Jtilde+bet*eye(Nel+1))\(Jtilde'*InvGamma_n*(Uel-Uref.Electrode(:))-bet*(beta-beta_0));
    
    
    
     deltatheta=(Jtilde'*InvGamma_n*Jtilde+IGamma_sigma_z)\(Jtilde'*InvGamma_n*(Uel-Uref.Electrode(:))-IGamma_sigma_z*(beta-beta_0))
    
    
    %     incl=find(deltatheta< -10);
    %      deltatheta(incl)=-10;
    %
    %          incl2=find(deltatheta> 1);
    %      deltatheta(incl2)= 1;
    % lograritminen
    
    % Log
    
    Uhat=Uref.Electrode;
    figure(3),clf
    plot(Uel(:,1),'b'),hold on
    plot(Uhat(:),'r')
    legend('Meas','Est')
    drawnow
    
    beta=beta+k*deltatheta;
    incl2=find(beta < -10);
    beta(incl2)= -10;
    incl3=find(beta(2:end) > 1);
    beta(incl3+1)= 1;
    
    z=exp(beta(2:(Nel+1)))
    sigma=exp(beta(1))
    z_ii(:,ii)=z(:);
    sigma_ii(ii)=sigma;
    
    
    if (ii<2)
        kriteeri(ii)=1;
    else
        kriteeri(ii) = norm(z_ii(:,ii)-z_ii(:,ii-1))
    end
    ii=ii+1;
end

% clear thta

index=find(kriteeri==min(kriteeri));

sigma_z(1)=sigma_ii(index);
sigma_z(2)=mean(z_ii(:,index));
end

