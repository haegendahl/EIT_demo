element_width = 1.1;
% assumed maximum change of conductivity within an elemement
sigdiffmamax = 0.5;
% approximate maximum of gradient
siggradmamax = sigdiffmamax/element_width;
% porcentage of confiability
p = 99.9;
% Optimal weigthing parameter
AKU_alpha = -log(1-p/100)/siggradmamax;


tic
%CASES = [1e3 1e2 AKU_alpha 1e-3 1e-4 1e-5 1e-6]
MakePaths
CASES = [1e-5]
RandStream.setDefaultStream(RandStream('mt19937ar','seed',10));

for gg=1:length(CASES)
    
%-------

clc
%clear all
clearvars -except CASES gg
close all



%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%

% - FORWARD MESH for Simulating the EIT Measurements

GEOMETRYTYPE = 'watertank3D_2'
%geofilename = 'watertank3D_2.geo'; % sphere inclusion
%geofilename = 'watertank3D_3.geo'; % cylinder inclusion
geofilename = 'watertank3D.geo'; % cylinder inclusion

Nel = 64;


FORWARDMESHREADY = 0;
fwd1stmeshname = 'TANK_MESH_1st';
fwd2ndmeshname = 'TANK_MESH_2nd';

order = 2;  


ng_parameters = {'Mesh granularity', 5,...   % overall mesh density, 1 = coarsest, 5 = finest
                 '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
                 'Max mesh-size', 2,...   % maximum allowed element size
                 'Mesh-size grading', .2,... 
                 'Elements per curvature', 2,... 
                 'Elements per edge', 2.0};
MaxElemSizeUnderElectrode = .5;

% - FORWARD MESH for solving the INVERSE PROBLEM (INVERSE MESH)

GEOMETRYTYPE_INV = 'watertank3D'
geofilename_inv = 'watertank3D.geo';

INVERSEMESHREADY = 0;
inv1stmeshname = 'TANK_MESH_1st_inv';
inv2ndmeshname = 'TANK_MESH_2nd_inv';

agradfilename = 'watertank3DAgrad';
order = 2;

ng_parameters_inv = {'Mesh granularity', 3,...   % overall mesh density, 1 = coarsest, 5 = finest
                 '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
                 'Max mesh-size', 2,...   % maximum allowed element size
                 'Mesh-size grading', .2,... 
                 'Elements per curvature', 2,... 
                 'Elements per edge', 2.0};
MaxElemSizeUnderElectrode = .5;


UNIFORMMESH = 0;
INTMATRIX_INV_1ST_READY = 0;
invmeshtype = 'Uniform';
elementsize_inv = 1.1;
invmeshfile = 'TANK_UNIFORM_INV';
intpmatfile_inv_1st = 'Mapping_to_unimesh_1st';



%%%%%%%%%%%%%%%%%%
% Forward meshes %
%%%%%%%%%%%%%%%%%%

if ~FORWARDMESHREADY

  [HF,gF,h2nd,NodeF,ElementF,elind2ndF,eltetra2ndF,E2ndF] = makecylmesh(GEOMETRYTYPE, geofilename,2,ng_parameters,fwd2ndmeshname,MaxElemSizeUnderElectrode,Nel);
  [HF1st,gF1st,NodeF1st,ElementF1st,elindF1st,eltetraF1st,EF1st] = Reduce2ndOrderMesh(HF,gF,elind2ndF,fwd1stmeshname);
  figure(2),simp_plot_3d(gF1st,HF1st), view(40,12)
  
  %% For the Jacobian
%   WW = speye(size(gF1st,1));
%   Agrad = jacobian3dnodeFast(Node1st,Element1st,WW);
%   eval(['save ',agradfilename,' Agrad'])

else

  mS = load(fwd1stmeshname);
  gF1st = mS.g;
  HF1st = mS.H;
  NodeF1st = mS.Node;
  ElementF1st = mS.Element;
  eltetraF1st = mS.eltetra;
  EF1st = mS.E;
  clear mS


  mS = load(fwd2ndmeshname);
  gF = mS.g;
  HF = mS.H;
  NodeF = mS.Node;
  ElementF = mS.Element;
  clear mS

end

% 2nd order organizing...
[ElementF] = organize_faces2nd(ElementF,gF); 
[HF ElementF] = organize_tetras2nd(HF,ElementF,gF);


%%%%%%%%%%%%%%%%%%
% Inverse mesh   %
%%%%%%%%%%%%%%%%%%

if ~INVERSEMESHREADY

  [H,g,h2nd,Node,Element,elind2nd,eltetra2nd,E2nd] = makecylmesh(GEOMETRYTYPE_INV, geofilename_inv,2,ng_parameters_inv,inv2ndmeshname,MaxElemSizeUnderElectrode,Nel);
  [H1st,g1st,Node1st,Element1st,elind1st,eltetra1st,E1st] = Reduce2ndOrderMesh(H,g,elind2nd,inv1stmeshname);
  figure(2),simp_plot_3d(g1st,H1st), view(40,12)
  
  %% For the Jacobian
  WW = speye(size(g1st,1));
  Agrad = jacobian3dnodeFast(Node1st,Element1st,WW);
  eval(['save ',agradfilename,' Agrad'])

else

  mS = load(inv1stmeshname);
  g1st = mS.g;
  H1st = mS.H;
  Node1st = mS.Node;
  Element1st = mS.Element;
  eltetra1st = mS.eltetra;
  E1st = mS.E;
  clear mS
  eval(['load ',agradfilename])

  mS = load(inv2ndmeshname);
  g = mS.g;
  H = mS.H;
  Node = mS.Node;
  Element = mS.Element;
  clear mS

end

% 2nd order organizing...
[Element] = organize_faces2nd(Element,g); 
[H Element] = organize_tetras2nd(H,Element,g);


%%%%%%%%%%%%%%%%%
% UNIFORM MESH  %
%%%%%%%%%%%%%%%%%

if ~UNIFORMMESH

  if strcmp(GEOMETRYTYPE_INV,'watertank3D'); % || strcmp(GEOMETRYTYPE_INV,'watertank3D_2')  ;
    cylinder_height = max(g(:,3));
    cylinder_radius = max(g(:,1));
    radius_div = round(cylinder_radius/elementsize_inv);
    
  end
  
  if(strcmp(invmeshtype,'Uniform'))
    [ginv,Hinv] = MakeUniformInvGrid(cylinder_height,cylinder_radius,radius_div,'dontdraw',invmeshfile);
  elseif(strcmp(invmeshtype,'Forwardmesh'))
    ginv = g1st;
    Hinv = H1st;
    P1st = eye(size(ginv,1));
  end
else
  eval(['load ',invmeshfile]) 
end
figure(3),clf,simp_plot_3d(ginv,Hinv)

Node_inv = [];
Nodeinv_P = MakeNode3dSmallFast(Hinv,ginv);
Element_inv = [];

if INTMATRIX_INV_1ST_READY
   if(~strcmp(invmeshtype,'Forwardmesh'))
      eval(['load ',intpmatfile_inv_1st])
   end
else
   sigmatmp = randn(size(ginv,1),1);
   [sigmatmp1st,P1st] = Interpolate2Newmesh3DNode(ginv,uint32(Hinv),Nodeinv_P,sigmatmp,g1st,[]); 
   eval(['save ',intpmatfile_inv_1st,' P1st']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATING THE TARGET %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_sigma = 0.5;
max_sigma = 0.5;



%theta_true = min_sigma*ones(size(gF1st,1),1); 

%Volumen Triangular

% theta_true = max_sigma*ones(size(gF1st,1),1); 
% 
% n1=0;
% n2=1;
% n3=0;
% 
% index = ((n1*gF1st(:,1)-5) +(1*gF1st(:,2)-0) + (n3*gF1st(:,3)-0)>= 0 & (1*gF1st(:,1)+5) +(-1*gF1st(:,2)+5) + (n3*gF1st(:,3)-0)>= 0 & (1*gF1st(:,1)-5) +(n2*gF1st(:,2)-7) + (n3*gF1st(:,3)-0)<= 0);


% crack
% 
theta_true = max_sigma*ones(size(gF1st,1),1); 
width= 2;
R = 12;
index = (gF1st(:,1)<=R & gF1st(:,1)>=0 & gF1st(:,2)<= width & gF1st(:,2)>= 0);

% cylinder
% x0=5;
% y0=3;
% r=3;
% index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2 <= (r+1e-1)^2;
% 

% sphere
% x0=0;
% y0=5;
% z0=10; 
% r=3;
% 
% index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2+(gF1st(:,3)-z0).^2 <= (r+1e-1)^2;

theta_true(index) = min_sigma;


%----

figure(1)
set(gcf,'renderer','zbuffer')
[hp1,hp2] = simp_plot_3dvec(gF1st,HF1st,'p(:,3)<5',theta_true,theta_true);
view(40,12)
axis tight
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate EIT measurements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 64                                        % Number of electrodes
                       % Number of electrodes
z= 0.00001*ones(L,1);          % Contact Impedances

disp(['Using: Opposite current injections and measurements',10]);
I = toeplitz([1; zeros(L/2-1,1); -1; zeros(L/2-1,1)]',[1 zeros(1,L/2-1)]); % Opposite current injections
MeasPatt = [ones(1,L-1);-eye(L-1)];                                        % Opposite measurements

Uel = ForwardSolution3d2ndElectrode_fix(NodeF,ElementF,I,theta_true,z,MeasPatt,'real');
Uel = Uel(:);

figure(1)
plot(Uel(:))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Adding Noise to the Measurements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE:

% Measurement Noise Statistics

meas_noise_coef = 2e-4;  % (meas_noise_coef*(max(Uel)-min(Uel)))^2       
meas_noise_coef2 = 0; %  var(l) = (meas_noise_coef2*Uel(l))^2
% meas_noise_coef_e = meas_noise_coef; % estimated noise level 
% gamma_kesk = 0; epsilon1 = 0; % for approximation error

beta_k = (meas_noise_coef*(max(max(Uel))-min(min(Uel))))^2;
var_Uel = beta_k + (meas_noise_coef2*abs(Uel)).^2;
gamma = diag(var_Uel(:));

% ---- Computing the Gamma^-1 (per block)

invGamma = sparse(inv(gamma));
Ln = chol(invGamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% homogenous estimate for the admittivity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse
ng1 = size(ginv,1);
nH1 = size(Hinv,1);
ng = ng1;

% forward
ng2 = size(g,1);
nH2 = size(H,1);

% 1st order (jacobian)
ng1st = size(g1st,1);
nH1st = size(H1st,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         INVERSE PROBLEM             %%%
%%%    Solving the Inverse Problem      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxIter= 25;
theta_err = 0;
%thetaexp = theta0;

% Estabilization parameter
beta    = 1e-4;

% TV Weighting parameter
% R and Ai - Matrices
% Note: now A_i is assigned as A_i = 1, for all i

R = MakeGradient3DFast_fix(ginv,Hinv);
Ai = ones(nH1,1);


% PlotGradient(theta_true,R_true, gF,HF,2),
% axis image,colorbar,
% axis off
% colorbar
% colormap(hot)
% title('Gradient')



%--- Selection of alpha  ----%


element_width = elementsize_inv;
% assumed maximum change of conductivity within an elemement
sigdiffmamax = max(theta_true) - min(theta_true);
% approximate maximum of gradient
siggradmamax = sigdiffmamax/element_width;
% porcentage of confiability
p = 99.9;
% Optimal weigthing parameter
alpha = -log(1-p/100)/siggradmamax;

alpha = CASES(gg)



%----  Prior for contact impedances ------%


zexpect = 10^(-1);         % expectation of z (if FixedExpectations)
z_bg_max_coef = 1.0001;
zmax_coef = 100;

zmax = zmax_coef*zexpect;
zbgmax = z_bg_max_coef*zexpect;


% contact impedance prior

z0 = 0;
zmin = z0 - zmax;
zrange = [zmin zmax];
Gamma_z = ((diff(zrange)/6)^2) * eye(Nel); 
% background

zbg0 = zexpect;
zbgmin = 2*zbg0 - zbgmax;
zbgrange = [zbgmin zbgmax];
bg_var = (diff(zbgrange)/6)^2;
Gamma_z = Gamma_z + bg_var*ones(Nel);
IGamma_z = inv(Gamma_z);
Lz =  chol(IGamma_z);

%%%-----------------------------

% Use this as initial guess in Gauss-Newton iteration

% first/zero iterate (initial guess)
theta0 = min_sigma; % - BACKGROUND
theta = theta0*ones(ng1,1);

z = zbg0*ones(Nel,1)

zest(:,1) = z; 

% ------ Starting the 0th iteration ------ %

disp(['Starting 0th iteration...',10])

admest = P1st*theta(:,end); 

Uref2 = ForwardSolution3dCI(Node1st,Element1st,I,admest,zest,MeasPatt,'real');
Current = Uref2.Current;
MsField = Uref2.MeasField;


zhat{1} = Uref2.Electrode;
zhat{2} = Current;
zhat{3} = MsField;

Jz = ContactImpedanceJacobi(g1st,H1st,z,zhat,eltetra1st,E1st);

% Jacobian

[ar ac av] = regroupAgrad2cell(Agrad,ng1st);

J = jacobian3dnodeCellReal(ar,ac,av,1,Uref2.Current(1:ng1st,:),Uref2.MeasField(1:ng1st,:),'real');
Jt = [J*P1st Jz];

Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,admest,zest,MeasPatt,'real');
Urefel = Urefel(:);

Fnorm(1) = 0.5 *norm(Ln*(Uel-Urefel))^2 + norm(Lz*(zest(:,1) - z))^2 + alpha*TV3D_functional(R, Ai, theta(:,end), beta);
measnorm(1) = norm(Uel-Urefel);
F0 = Fnorm(1);


disp(['Computing the Gradient of TV', 10]);
%grad_TV = Gradient_3DTV(R, Ai, theta(:,end), beta);
grad_TV = Gradient_3DTV_AKU(R, theta(:,end), beta);

disp(['Computing the Hessian of TV', 10]);
Hess_TV = Hessian_3DTV(R, Ai, theta(:,end), beta);

FIRSTSTEP = 0.5; % first step in line search
NMAX = 5;

disp(['Starting Gauss-Newton iterations...',10])


Wn = Ln'*Ln;
Wz = Lz'*Lz;



for ii = 2:maxIter
    
%     HH = (J'*W*J + alpha* Hess_TV);
%     zz =  J'*W* (Uel-Urefel) - alpha*grad_TV;
    
    HH =  Jt'*Wn*Jt + blkdiag(alpha*Hess_TV, 2*Wz);  
    zz =  Jt'*Wn*(Uel-Urefel) +[-alpha*grad_TV; 2*Wz'*(z-zest(:,ii-1))];
    
    dtheta=HH\zz;
    
    disp(['Computing k using Linesearch... ',10])
    F0 = Fnorm(ii-1);
    
    %   Linesearch usa bases en 2do orden para calcular el Potencial
    [k, ~, ~] = linesearch2(FIRSTSTEP,NMAX, dtheta, Uel,F0, Node, Element,I,theta(:,ii-1), P1st, z, MeasPatt, Ln, R, Ai, beta,alpha, 501);
    %k = 0.5
    if(k<=0)
        break;
    end
    
    disp(['STEP: ' num2str(k)])
    
    % new estimate
    theta(:,ii) = theta(:,ii-1) + k*dtheta;
    theta(theta(:,ii)<0,ii) = 1e-8;
        
    admest = theta(:,end);
        
    disp(['Computing the Gradient of TV', 10]);
    %grad_TV = Gradient_3DTV(R, Ai, theta(:,ii), beta);
    grad_TV = Gradient_3DTV_AKU(R, Ai, theta(:,end), beta);
    
    disp(['Computing the Hessian of TV', 10]);
    Hess_TV = Hessian_3DTV(R, Ai, theta(:,ii), beta);
    
    Urefel = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*admest,z,MeasPatt,'real');
    Urefel = Urefel(:);
    
    % Compute the Jacobian
    Uref2 = ForwardSolution3dCI(Node1st,Element1st,I,P1st*admest,z,MeasPatt,'real');
    J = jacobian3dnodeCellReal(ar,ac,av,1,Uref2.Current(1:ng1st,:),Uref2.MeasField(1:ng1st,:),'real');
    J = J*P1st;
    
    Fnorm(ii) = 0.5 *norm(Ln*(Uel-Urefel))^2  + alpha*TV3D_functional(R, Ai, theta(:,end), beta);
    measnorm(ii) = norm(Uel-Urefel)^2;
    

    ii

    theta_err(ii) = norm(theta(:,ii)-theta(:,ii-1))/max([1 norm(theta(:,ii))]);
    
    
    figure(5),clf,plot(Uel,'b'),hold on,plot(Urefel,'r')
    legend('Uel','Urefel')
    
    if(theta_err(ii) < 1e-3)
        disp('Break!')
        break
    end
end




%%%%%%%%%%%%%%%%%%%%%%
%%% Graphics Output %%
%%%
%%% Visualization  %%%
%%%%%%%%%%%%%%%%%%%%%%

IntpMatrices = [];
xyslices = max(ginv(:,3))*[0:.1:1]; 
xzslices = max(ginv(:,2))*[];
yzslices = max(ginv(:,1))*[];


cmax = max([theta_true(:); theta(:,end)]);
cmin = min([theta_true(:); theta(:,end)]);

range= [min_sigma, max_sigma];


[VisElementsF, VisNodesF, IntpMatricesF] = Generate2DSlices(gF1st,HF1st,NodeF1st,theta_true,xyslices,xzslices,yzslices,max(gF(:,1))/50,IntpMatrices);
[VisElements, VisNodes, IntpMatrices] = Generate2DSlices(ginv,Hinv,Nodeinv_P,theta(:,end),xyslices,xzslices,yzslices,max(ginv(:,1))/50,IntpMatrices);

PlotMyResults(VisElementsF, VisNodesF, IntpMatricesF, theta_true,range,101)
PlotMyResults(VisElements, VisNodes, IntpMatrices, theta(:,end),range,101)



time = clock;
foldername = [num2str(time(1)) '-' num2str(time(2)) '-' num2str(time(3)) '-' num2str(time(4)) '-' num2str(time(5)) '-' num2str(round(time(6)))];

mkdir(foldername)


%save figures here
PrintMyResults(VisElementsF, VisNodesF,IntpMatricesF,theta_true, 'XYslices_true',range, ['./',foldername]);
printcolorbar(min_sigma, max_sigma,['./',foldername,'/colorbar']);


PrintMyResults(VisElements, VisNodes,IntpMatrices,theta(:,end), 'XYslices',range, ['./',foldername]);
printcolorbar(min_sigma, max_sigma, ['./',foldername,'/colorbar']);

save(['./',foldername,'/data.mat'])
copyfile('./MAIN2.m',['./',foldername])



end
tic



% 
% 
% 
% 
% 
% viewangle = [-47 -10];
% grid on
% cmin = min_sigma;
% cmax = max_sigma;
% 
% 
% DrawTarget(g1st,H1st,Node1st,theta_true,101,102, viewangle,cmin,cmax);
% 
% DrawTarget(ginv,Hinv,Nodeinv_P,theta(:,end),201,202, viewangle,cmin,cmax);
% DrawCircle(0,0,0,12.75);
% 
% line([8.04,8.04],[-9.896,-9.896 ],[0,20.4],'Color','black','LineWidth',1);
% line([-8.506,-8.506],[9.498,9.498 ],[0,20.4],'Color','black','LineWidth',1);
% 
% DrawCircle(0,0,20.4,12.75);
% 
% 
% 
