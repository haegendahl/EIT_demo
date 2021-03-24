
% Impletation of EIT solver using 3D-TV prior. 
% Barrier function was used in the optimization problem.
% Based on original codes of M. Vauhkonen and K. Karhunen 
% 2013, G. Gonzalez

% Requires NETGEN 4.9.X for generating the meshes.

clc
clear all
close all

MakePaths

%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%

% - FORWARD MESH for Simulating the EIT Measurements

GEOMETRYTYPE = 'watertank3D_2'
%geofilename = 'watertank3D_2.geo'; % sphere inclusion
geofilename = 'watertank3D_3.geo'; % cylinder inclusion
%geofilename = 'wt_triangular_prism.geo'; % triangular prism


Nel = 64;

FORWARDMESHREADY = 1;
fwd1stmeshname = 'TANK_MESH_1st';
fwd2ndmeshname = 'TANK_MESH_2nd';

order = 2;

ng_parameters = {'Mesh granularity', 5,...   % overall mesh density, 1 = coarsest, 5 = finest
    '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
    'Max mesh-size', 0.7,...   % maximum allowed element size
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

ng_parameters_inv = {'Mesh granularity', 4,...   % overall mesh density, 1 = coarsest, 5 = finest
    '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
    'Max mesh-size', .8,...   % maximum allowed element size
    'Mesh-size grading', .2,...
    'Elements per curvature', 2,...
    'Elements per edge', 2.0};
MaxElemSizeUnderElectrode = .5;


UNIFORMMESH = 0;
INTMATRIX_INV_1ST_READY = 0;
invmeshtype = 'Uniform';
elementsize_inv = 1.5;
invmeshfile = 'TANK_UNIFORM_INV';
intpmatfile_inv_1st = 'Mapping_to_unimesh_1st';
intpmatfile_fwd_1st = 'Mapping_fwd_to_unimesh_1st';


%%%%%%%%%%%%%%%%%%
% Forward meshes %
%%%%%%%%%%%%%%%%%%

if ~FORWARDMESHREADY
    
    [HF,gF,h2ndF,NodeF,ElementF,elind2ndF,eltetra2ndF,E2ndF] = makecylmesh(GEOMETRYTYPE, geofilename,2,ng_parameters,fwd2ndmeshname,MaxElemSizeUnderElectrode,Nel);
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
        eval(['load ',intpmatfile_fwd_1st])
    end
else
    sigmatmp = randn(size(ginv,1),1);
    sigmatmp2 = randn(size(gF1st,1),1);
    [sigmatmp1st,P1st] = Interpolate2Newmesh3DNode(ginv,uint32(Hinv),Nodeinv_P,sigmatmp,g1st,[]);
    [sigmaFtmp1st,PF1st] = Interpolate2Newmesh3DNode(gF1st,uint32(HF1st),NodeF1st,sigmatmp2,ginv,[]);
    eval(['save ',intpmatfile_inv_1st,' P1st']);
    eval(['save ',intpmatfile_fwd_1st,' PF1st']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATING THE TARGET %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_sigma = 0.5;
max_sigma = 2;




%Volumen Triangular

%index = (0*gF1st(:,1)-0) +(-1*gF1st(:,2)-0) + (0*gF1st(:,3)-0)<= 1e-2 &  (-1*gF1st(:,1)-0) +(0*gF1st(:,2)-0) + (0*gF1st(:,3)-0)<= 1e-2 & (1*gF1st(:,1)-10) +(1*gF1st(:,2)-0) + (0*gF1st(:,3)-0)<= 1e-2;



% crack
%
% theta_true = max_sigma*ones(size(gF1st,1),1);
% width= 2;
% R = 12;
% index = (gF1st(:,1)<=R & gF1st(:,1)>=0 & gF1st(:,2)<= width & gF1st(:,2)>= 0);

%cylinder
x0=5;
y0=3;
r=3;
index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2 <= (r+1e-2)^2;


% sphere
% x0=0;
% y0=5;
% z0=10;
% r=3;
%
% index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2+(gF1st(:,3)-z0).^2 <= (r+1e-1)^2;

%theta_true = min_sigma*ones(size(gF1st,1),1);
theta_true = max_sigma*ones(size(gF1st,1),1);
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
z= 0.00001*ones(L,1);          % Contact Impedances: This value is a strong assumption. 

disp(['Using: Custom current injections and adj measurement pattern',10]);
I1 = -eye(L/4,L/4);

I2 = zeros(L/2,L/4);

I3 = toeplitz([zeros(1,(L/8)), 1  zeros(1,L/8-1)]);

I = [I1;I2;I3]; % First layer oppsite last layer

%MeasPatt = [ones(1,L-1);-eye(L-1)];                                        % Opposite measurements

MeasPatt = Current(64,0,'adj');                                       % Opposite measurements

U = ForwardSolution3d2ndElectrode_fix(NodeF,ElementF,I,theta_true,z,MeasPatt,'real');
Uel_noiseless = U(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Adding Noise to the Measurements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE:

% Measurement Noise Statistics

meas_noise_coef = 2e-4;
meas_noise_coef2 = 1e-3;

beta_k = (meas_noise_coef*(max(max(Uel_noiseless))-min(min(Uel_noiseless))))^2;
var_Uel = beta_k + (meas_noise_coef2*abs(Uel_noiseless)).^2;
gamma_n = diag(var_Uel(:));

invGamma = sparse(inv(gamma_n));
Ln = chol(invGamma);

% add noise to the measurements

Uel1 = Uel_noiseless + meas_noise_coef2*abs(Uel_noiseless).*randn(size(Uel_noiseless));
Uel = Uel1 + (meas_noise_coef*(max(Uel1)-min(Uel1)))*randn(size(Uel1));


figure(103)
plot(Uel_noiseless,'r')
hold on
plot(Uel)
legend('Noiseless','Noisy')


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


%%%%%%%%%%%%%%%%%%%%%%%
%%% INVERSE PROBLEM %%%
%%%%%%%%%%%%%%%%%%%%%%%

[ar ac av] = regroupAgrad2cell(Agrad,ng1st);
clear Agrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Solving the Inverse Problem      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxIter= 25;
theta_err = 0;
%thetaexp = theta0;

% Estabilization parameter
beta = 1e-3;

% Barrier method, penalty parameter
% changes adaptively
% Will be initialized later...
a = 0;

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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Selection of alpha %%
%%%%%%%%%%%%%%%%%%%%%%%%%

element_width = mean_edge(ginv,Hinv);
% assumed maximum change of conductivity within an elemement
sigdiffmamax = max(theta_true) - min(theta_true);
% approximate maximum of gradient
siggradmamax = sigdiffmamax/element_width;
% porcentage of confiability
p = 99;
% Optimal weigthing parameter
AKU_alpha = -log(1-p/100)/siggradmamax;

alpha = AKU_alpha;


% Use this as initial guess in Gauss-Newton iteration

% first/zero iterate (initial guess)
theta0 = min_sigma; % - BACKGROUND
theta = theta0*ones(ng1,1);


% ------ Starting the 0th iteration ------ %

disp(['Starting 0th iteration...',10])

admest = theta(:,end);

% - Interior Point Method Patch

FIRSTSTEP = 0.1; % first step in line search
maxIter = 25;
NMAX = 10;
minval = 10^(-6);
DRAW   = 1;


cost_params.Node = Node;
cost_params.Element =  Element;
cost_params.I = I;
cost_params.P1st = P1st;
cost_params.z = z;
cost_params.MeasPatt = MeasPatt;
cost_params.Uel = Uel;
cost_params.Ln = Ln;
cost_params.R = R;
cost_params.Ai = Ai;
cost_params.beta = beta;
cost_params.alpha = alpha;
cost_params.a = a;



% - Prior Model

Ffunc = 'costfunTVnoneg';
[Fnorm(1),Urefel,Uref] = costfunTVnoneg(theta(:,end),cost_params);
disp(['Computing the Gradient of TV', 10]);
grad_TV = Gradient_3DTV(R, theta(:,end), beta);
disp(['Computing the Hessian of TV', 10]);
Hess_TV = Hessian_TV3D_fast_2(R,theta(:,end),beta);


% Barrier method, initialization of penalty parameter a

%c = .0001;
c = 1e-6
l = sum(log(theta(:,end)));
q = -a*l;
Fnorm(1) = Fnorm(1) - q;
lm = -log(minval);
a = (c/lm)*Fnorm(1);
a_ALL = a
a=0;
cost_params.a = a;
q = -a*sum(log(theta(:,end)));
Fnorm(1) = Fnorm(1) + q;
cost_params.a = a;

gradq = -a./theta(:,end);
Hessq = diag(a./(theta(:,end)).^2);

measnorm(1) = norm(Uel-Urefel);
F0 = Fnorm(1);


% Jacobian
Uref2 = ForwardSolution3dCI(Node1st,Element1st,I,P1st*theta(:,end),z,MeasPatt,'real');
J = jacobian3dnodeCellReal(ar,ac,av,1,Uref2.Current(1:ng1st,:),Uref2.MeasField(1:ng1st,:),'real');
J = J*P1st;


disp(['Starting Gauss-Newton iterations...',10])

W = Ln'*Ln;

for ii = 2:maxIter
    
    HH = (J'*W*J + alpha* Hess_TV + Hessq);
    zz =  J'*W* (Uel-Urefel) - alpha*grad_TV - gradq;
    dtheta=HH\zz;
    
    disp(['Computing k using Linesearch... ',10])
    F0 = Fnorm(ii-1);
    
    %   Linesearch usa bases en 2do orden para calcular el Potencial
    [k, ~, ~] = linesearch_EIT(theta(:,ii-1),dtheta,F0,Ffunc,FIRSTSTEP,minval,NMAX,501,DRAW,cost_params);
    
    
    if(k<=0)
        break;
    end
    
    disp(['STEP: ' num2str(k)])
    
    % new estimate
    theta(:,ii) = theta(:,ii-1) + k*dtheta;
    theta(theta(:,ii)<0,ii) = minval;
    
    admest = theta(:,end);
    

    [Fnorm(ii),Urefel,Uref] = costfunTVnoneg(theta(:,end),cost_params);

    disp(['Computing the Gradient of TV', 10]);
    grad_TV = Gradient_3DTV(R, theta(:,ii), beta);
    disp(['Computing the Hessian of TV', 10]);
    Hess_TV = Hessian_TV3D_fast_2(R,theta(:,ii), beta);

    
    
    % Barrier method
    % reselect the penalty parameter a
    
    l = sum(log(theta(:,ii)));
    q = -a*l;
    Fnorm(ii) = Fnorm(ii) - q;
    a = (c/lm)*Fnorm(ii);
    a_ALL = [a_ALL a]
    cost_params.a = a;
    % a = 0;
    q = -a*sum(log(theta(:,ii)));
    Fnorm(ii) = Fnorm(ii) + q;
    
    gradq = -a./theta(:,ii);
    Hessq = diag(a./(theta(:,ii)).^2);
    
    measnorm(ii) = norm(Uel-Urefel)^2;
    
    
    
    % Compute the Jacobian
    Uref2 = ForwardSolution3dCI(Node1st,Element1st,I,P1st*theta(:,end),z,MeasPatt,'real');
    J = jacobian3dnodeCellReal(ar,ac,av,1,Uref2.Current(1:ng1st,:),Uref2.MeasField(1:ng1st,:),'real');
    J = J*P1st;
    
    
    figure(1),clf,plot(Uel,'b'),hold on,plot(Urefel,'r')
    legend('Uel','Urefel')
    
    
    figure(2), clf
    set(gcf,'renderer','zbuffer')
    subplot(1,2,1)
    [hp1,hp2] = simp_plot_3dvec(ginv,Hinv,'p(:,3)<=8',admest,admest);
    colormap(fireprint());
    colorbar
    %axis tight, view(-40,60)
    axis tight, view(0,90)
    subplot(1,2,2)
    [hp1,hp2] = simp_plot_3dvec(ginv,Hinv,'p(:,2)>=0',admest,admest);
    colormap(fireprint());
    axis tight, view(3)
    drawnow

    
    ii
    
    theta_err(ii) = norm(theta(:,ii)-theta(:,ii-1))/max([1 norm(theta(:,ii))]);

    
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
close all

% Generating folder tag.

time = clock;
RUN_TAG = [num2str(time(1)) '-' num2str(time(2)) '-' num2str(time(3)) '-' num2str(time(4)) '-' num2str(time(5)) '-' num2str(round(time(6)))];


MESH_COUNTER   = hh;
BCASES_COUNTER = ff;
CASES_COUNTER  = gg;

foldername = ['MESH', num2str(hh)]

%mkdir(foldername)

subfoldername = ['beta', num2str(ff)]

subsubfoldername =  ['c', num2str(gg)]

pathtarget = ['./',foldername,'/', subfoldername,'/', subsubfoldername];

mkdir(pathtarget);

type_of_slice_F = 'F_XYslices';
type_of_slice = 'XYslices';

printtype = 'png';


thta = PF1st*theta_true;

slicesF = GetSlices(ginv,thta);
slices  = GetSlices(ginv,theta);


% Graphic output

% - True theta


PlotScaledSlices(slicesF, pathtarget,type_of_slice_F,'Yes');

printcolorbar(min_sigma, max_sigma, [pathtarget,'/colorbar']);

PlotGradientSlices(slicesF, pathtarget,type_of_slice_F,'Yes');


% - Estimated theta

PlotScaledSlices(slices, pathtarget,type_of_slice,'Yes');

printcolorbar(min_sigma, max_sigma, ['./',foldername,'/colorbar']);

PlotGradientSlices(slices, pathtarget,type_of_slice,'Yes');



copyfile('./autocrop.*',[pathtarget])
%system([pathtarget,'/autocrop.py'])


save([pathtarget ,'/data.mat'])
copyfile('./MAIN_INT_POINT_METH.m',[pathtarget])


fid = fopen([pathtarget ,'/alpha=',num2str(CASES(gg))],'w');
fclose(fid)

fid = fopen([pathtarget ,'/beta=',num2str(BCASES(ff))],'w');
fclose(fid)

fid = fopen([pathtarget ,'/element_width=',num2str(element_width)],'w');
fclose(fid)

fid = fopen([pathtarget ,'/',num2str(RUN_TAG)],'w');
fclose(fid)




