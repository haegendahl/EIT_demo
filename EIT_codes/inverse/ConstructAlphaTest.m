function [eta,iGp] = ConstructAlphaTest(g,H,Node,Element,cond,alpha,range,indexes,PriorStruct)
%keyboard;
alpha0 = alpha;

%disp('Plane anisotropic prior!')
%R = regmatnode2D(g,H);

disp(['Using: ' PriorStruct.priortype])
eval(['R =' PriorStruct.priortype '(Node,Element,PriorStruct);']);
%R = MakeGradientRegMat3DFast(Node,Element);
% R = MakeRadialAnisotropicRegMat3D(Node,Element,[1 10 10],7.5);
%R = MakeAnisotropicRegMat3D_kiekot(Node,Element,[1 3 5],7.5);
%R = MakeAnisotropicRegMat3D_kiekot_horizontalbar(Node,Element,[1 5 3],7.5);

%R = MakeAnisotropicRegMat3DPlane(Node,Element,[0 50;0 50;0 4],[50 50 50]);
%R = MakeAnisotropicRegMat3DPlaneCrack(Node,Element,[0 50;0 50;0 4],[50 50 50]);
%R = alpha0*R;
cond0 =cond;

sN = size(g,1);
ind = 1:sN; 
%% index set I2

%ind1 = find( (g(:,1) == xmax |  g(:,2) == ymin) | (g(:,1) == xmax | g(:,2) == ymax));
%ind11 = find( (g(:,2) == xmin |  g(:,1) == ymin) | (g(:,2) == xmin | g(:,1) == ymax) );
%ind1 = [ind1;ind11];
%ind1 = unique(ind1); % indices of f1 (conditioned pixels)

%ind1 = ind(1:60:end);
ind1 = indexes;
ind2 = setdiff(ind,ind1); % indices of f2

%iG2 = alpha0*(R'*R);
iG2 = (R'*R);
  
G22 = iG2(ind2,ind2);
G12 = iG2(ind1,ind2);
G21 = iG2(ind2,ind1);
  
%P = (1/(.015*cond0)^2) * eye(length(ind1));
s0 = (diff(range)/6)^2

P = (1/s0) * eye(length(ind1));
iP = s0*eye(length(ind1));

% Construct the joint covariance. Note the
% correct replacement of the matrix blocks
iG22 = inv(G22);
iGp = zeros(sN,sN);

GG = iG22*G21;
d1 = diag(iG22);
d2 = diag(GG*iP*GG');
nd = length(d1);

% mean alpha criterion
m1 = sum(d1);
m2 = sum(d2);
alpha0 = (m1/(2*nd*s0 - 2*m2))
%alpha0 = 20000

% half maximum alpha criterion (nonlinear, bad)
%m1 = max(d1);
%m2 = max(d2);
%alpha0 = (m1/(4*s0-2*m2))

iGp(ind1,ind1) = 2*alpha0*G12*GG+P;
iGp(ind1,ind2) = 2*alpha0*G12;
iGp(ind2,ind1) = 2*alpha0*G21;
iGp(ind2,ind2) = 2*alpha0*G22;

% Construct the midpoint
eta = zeros(sN,1);
eta(ind1) = cond0;
eta(ind2) = -GG*eta(ind1);
 


%trimesh(H,g(:,1),g(:,2))
%hold on
%plot(g(ind1,1),g(ind1,2),'o')
