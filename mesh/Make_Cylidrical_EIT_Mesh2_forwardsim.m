% Aku Seppänen 7.2.2003
%
% Constructs two grids for cylindrical pipe
% 1) g,H: inverse computations
% 2) g2,H2: Forward simulations
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for inverse problem: Convection-diffusion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

el_width = 1.3;
bel_width = .76;
el_width_z = 1.3;
L = 16;

%zs = el_width_z/2*ones(2*(2*4+4+3*2),1);
%r = L*(el_width+bel_width)/(2*pi)*[1 .87 .69 .52 .36  0.18,  0];
%NN = [48 27 24 21 14 6 1];

zs2 = el_width_z/2*ones(2*(2*4+4+4*2),1);
r2 = L*(el_width+bel_width)/(2*pi)*[1 .85 .67 .5 .34  0.15,  0];
NN2 = [48 24 20 22 16 8 1];



%el = [2,1]; 
%h = sum(zs);
%N = [0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0  0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0];

el2 = [2,1]; 
h2 = sum(zs2);
N2 = [0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0];

[Eind2,elind2,g2,H2] = MakeMesh3D_FlowLoop(N2,r2,NN2,h2,el2,zs2,L,el_width,bel_width);

figure(1),clf
trimesh(H2,g2(:,1),g2(:,2),ones(size(g2,1),1),ones(size(g2,1),1)),axis image
view(2)

sN2 = max(size(g2));

%figno = 1;
%Draw_Cylindrical_Grid_new(g,H,elind,r,figno)

figno = 2
Draw_Cylindrical_Grid_new(g2,H2,elind2,r2,figno)


[Element2]=MakeElement3dSmall(H2,elind2,Eind2);
[Node2]=MakeNode3dSmall(Element2,g2);

sN2 = max(size(g2));
WW2 = speye(sN2,sN2);
Agrad2 = jacobian3dnode(Node2,Element2,WW2);


eval(['save ',meshfile_forwardsim,' Node2 Element2 g2 H2 N2 elind2 Eind2 Agrad2 zs2 NN2 sN2'])       

