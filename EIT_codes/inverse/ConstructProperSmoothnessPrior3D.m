function [eta,iGp] = ConstructProperSmoothnessPrior3D(g,H,Node,Element,cond,alpha,range,indexes)
%keyboard;
alpha0 = alpha;

%disp('Plane anisotropic prior!')
%R = regmatnode2D(g,H);

R = MakeGradientRegMat3DFast(Node,Element);
% R = MakeRadialAnisotropicRegMat3D(Node,Element,[1 10 10],7.5);
%R = MakeAnisotropicRegMat3D(Node,Element,[1 50 50],7.5);
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
P = (1/(diff(range)/6)^2) * eye(length(ind1));

% Construct the joint covariance. Note the
% correct replacement of the matrix blocks
iG22 = inv(G22);
iGp = zeros(sN,sN);
    
iGp(ind1,ind1) = 2*alpha0*G12*iG22*G21+P;
iGp(ind1,ind2) = 2*alpha0*G12;
iGp(ind2,ind1) = 2*alpha0*G21;
iGp(ind2,ind2) = 2*alpha0*G22;

% Construct the midpoint
eta = zeros(sN,1);
eta(ind1) = cond0;
eta(ind2) = -iG22*G21*eta(ind1);
 


%trimesh(H,g(:,1),g(:,2))
%hold on
%plot(g(ind1,1),g(ind1,2),'o')
