% Aku Seppänen 5.4.2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MakePaths




     %%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%% Grid %%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%


%%% Mesh for inverse problems %%%

zs=1.0*ones(11,1);
r=5.5*[1,0.8,0.55,0.3,0];
NN=[32 32 16 8 1];
el=[1,1]; 
h = sum(zs);
N = [0 0 1 0 0 1 0 0 1 0 0];

[Eind,elind,g,H] = MakeMesh3D(N,r,NN,h,el,zs);

figure,
trimesh(H,g(:,1),g(:,2),ones(size(g,1),1),ones(size(g,1),1)),axis image
view(2)

R = max(r);
Draw_Cylindrical_Grid

L = 3*16;                           % Number of electrodes
z = 0.01*ones(L,1);              % contact impedances
sN = max(size(g));             % number of parameters 
WW = speye(sN,sN);

[Element]=MakeElement3dSmall(H,elind,Eind);
[Node]=MakeNode3dSmall(Element,g);

Agrad = jacobian3dnode(Node,Element,WW);





%%% Grid for 2nd order basis %%%

[g4,H4,elind4,Eind4] = mesh2nd(g,H,elind,Eind,L);


%%%%%%% hilan tulostus %%%%%%%%%%%%%%%
clf
plot3(g4(:,1),g4(:,2),g4(:,3),'w.');hold on
nH=length(H4) %;koe=[];

figure(1)
hold on
for ii=1:length(H4)
  Hii=g4(H4(ii,:),:);
  Hii=[Hii(1,:);Hii(5,:);Hii(2,:);Hii(6,:);Hii(3,:); ...
     Hii(7,:);Hii(1,:);Hii(8,:);Hii(4,:); ...
     Hii(9,:);Hii(2,:);Hii(6,:);Hii(3,:);Hii(10,:);Hii(4,:)];
hHii=plot3(Hii(:,1),Hii(:,2),Hii(:,3),'k');
%drawnow
end

view(3)


%%%%%%%%%%%% elektrodien tulostus %%%%%%%%%%%%%%
for i=1:length(elind4)
 %t=[1 4 2 5 3 6];
 t=[1 2 3 7 8 9];
 u=fill3(g4(elind4(i,t),1),g4(elind4(i,t),2),g4(elind4(i,t),3),'r');
end

 view(3)

%%% Tallentaminen struktuureina

elind4b=elind4;
Eind4b=Eind4;
for jj=1:L   
 eval(['elind4.ElecNo', int2str(jj),'=elind4b(jj,:);']) 
 eval(['Eind4.ElecNo', int2str(jj),'=Eind4b(jj,:);'])
end



     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%% Conductivity distribution %%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     sigma_min = .5*1/200;        % Minimum value of sigma
     sigma_max = 1/200;           % Maximum value of sigma

     wdth = .3;
     rr = (max(g(:,1))-min(g(:,1)))/2;
     sigma = GenerateBubblesinGrid3D(g,sigma_min,sigma_max,wdth,4,rr);

     N_frames = 1;                % Number of frames of measurements
     N_c_inj = 16;                % Number of current injection patterns per frame
     TTT = N_frames*N_c_inj;      % total number of current injections



A = FEMmatrix3d2ndLinSigma(g4,H4,Eind4,elind4,z,sigma);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%% Parameters of the voltage measurements %%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     I_rms = 1;                   % rms-value of the electric current
     %meas_type = 'adj_meas';      % adjacent measurements
     meas_type = 'opp_meas';      % opposite measurements
     inj_type = 'across_pipe'
     %inj_type = 'acros_pipe2'
     inj_type = 'adjacentinj';

     Set_Meas_Parameters

     meas_noise_coef = .0001;      % Noise coefficient: variance of the noise is
                                  %    (meas_noise_coef*(max(Uel)-min(Uel)))^2
                                  %    e.g. .01 corresponds to 'one percent noise level'
     meas_noise_coef2 = .01;     % Noise coefficient: variance of the noise is
     %meas_noise_coef2 = .00;      %    var(l) = (meas_noise_coef2*Uel(l))^2
                                  %    e.g. .01 corresponds to 'one percent noise level'

 
     sigma_pr = 1/200;           % initial quess for the average conductivity



[U,p,r]=ForwardSolution3dnode(Node,Element,T,sigma,z,MeasPatt,'real'); 
Uel_nonoise = U.Electrode(:);




[I,C]=current_eit(L,1,'cee',g4);
C = C(:,1:L-1);
II1=MeasPatt'*sparse([zeros(L,size(g4,1)),C]);
II1=II1';
II=[zeros(size(g4,1),size(T,2));C'*T];
p=symmmd(A);
r1(p)=1:max(size(p));
R=chol(A(p,p));     
UU=R\(R'\[II1(p,:),II(p,:)]);
UU=UU(r1,:);
Uel=MeasPatt'*C*UU(size(g4,1)+1:size(A,1),size(II1,2)+1:size(UU,2));
Uel = Uel(:);

plot(Uel(1:100))
hold on
plot(Uel_nonoise(1:100),'r')
