% Aku Seppänen 7.2.2003
%
% Constructs two grids for cylindrical pipe
% 1) g,H: inverse computations
% 2) g2,H2: Forward simulations
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for inverse problem: Convection-diffusion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zs4 = 1/1.5*ones(1.5*20,1);
r4 = 5.5*[1 .8 .6 .4  .2,  0];
NN4 = round(.9*2*pi*r4);  NN4(length(NN4)) = 1;

el4 = [2,2]; 
h4 = sum(zs4);
N4 = zeros(1,size(zs4,1));N4(1)=1;

%r=5.5*[1,0.85,0.6,0.3,0];
%NN=[64 32 16 8 1];

h4 = sum(zs4);
N4 = zeros(1,size(zs4,1));N4(1)=1;

[Eind4,elind4,g4,H4] = MakeMesh3D(N4,r4,NN4,h4,el4,zs4);

figure,
trimesh(H4,g4(:,1),g4(:,2),ones(size(g4,1),1),ones(size(g4,1),1)),axis image
view(2)

sN4 = max(size(g4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for inverse problem: EIT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%zs=1.0*ones(20,1);
%r=5.5*[1,0.85,0.6,0.3,0];
%NN=[64 32 16 8 1];
%el=[2,2]; 
%h = sum(zs);
%N = [0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 0 0];

zs = 1/1.5*ones(1.5*20,1);
r = 5.5*[1 .92 .8 .6 .4  .2,  0];
NN = round(.9*2*pi*r); NN(length(NN)) = 1;
NN(1) = 64; NN(2) = 32;
el=[2,2]; 
h = sum(zs);
N = [0 0 0 0 0 1  1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0  0 1 1 0 0 0 0 0 ];


[Eind,elind,g,H] = MakeMesh3D(N,r,NN,h,el,zs);

figure,
trimesh(H,g(:,1),g(:,2),ones(size(g,1),1),ones(size(g,1),1)),axis image
view(2)

Draw_Cylindrical_Grid_Rotate

[Element]=MakeElement3dSmall(H,elind,Eind);
[Node]=MakeNode3dSmall(Element,g);

sN = max(size(g));
WW = speye(sN,sN);
Agrad = jacobian3dnode(Node,Element,WW);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for EIT simulations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%zs2 = 1.0*ones(20,1);
%r2 = 5.5*[1, 0.85, 0.65, 0.45, 0.25, 0];
%NN2 = [64 48 24 12 6 1];
%el2 = [2,2]; 
%h2 = sum(zs2);
%N2 = [0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 0 0];

zs2 = 1/1.5*ones(1.5*20,1);
r2 = 5.5*[1, 0.94, .8 0.65, 0.45, 0.25, 0];

r2 = 5.5*[1 .94 .85 .67 .5 .35  0.17,  0];
NN2 = round(1.2*2*pi*r2);  NN2(length(NN2)) = 1;
NN2(1) = 80; NN2(2) = 64;
%NN2 = [80 64 48 24 12 6 1];
el2 = [3,2]; 
h2 = sum(zs2);
N2 = [0 0 0 0 0 1  1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0  0 1 1 0 0 0 0 0 ];


[Eind2,elind2,g2,H2] = MakeMesh3D(N2,r2,NN2,h2,el2,zs2);

figure,
trimesh(H2,g2(:,1),g2(:,2),ones(size(g2,1),1),ones(size(g2,1),1)),axis image
view(2)

[Element2]=MakeElement3dSmall(H2,elind2,Eind2);
[Node2]=MakeNode3dSmall(Element2,g2);

sN2 = max(size(g2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for Convection-Diffusion computation: true evolution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zs3 = 1/1.5*ones(1.5*20,1);
%r3 = 5.5*[1, 0.85, 0.65, 0.45, 0.25,  0];
%NN3 = [40 36 24 12 6 1];

r3 = 5.5*[1 .85 .67 .5 .35  0.17,  0];
NN3 = round(1.2*2*pi*r3);  NN3(length(NN3)) = 1;

el3 = [2,2]; 
h3 = sum(zs3);
N3 = zeros(1,size(zs3,1));N3(1)=1;


[Eind3,elind3,g3,H3] = MakeMesh3D(N3,r3,NN3,h3,el3,zs3);

figure,
trimesh(H3,g3(:,1),g3(:,2),ones(size(g3,1),1),ones(size(g3,1),1)),axis image
%trimesh(H3,g3(:,1),g3(:,2),g3(:,3),ones(size(g3,1),1)),axis image
view(2)

sN3 = max(size(g3));





save /wrk3/3Ddyn/mesh/spipeGRID_NEW Node Element g H N elind Eind Agrad zs NN sN...
                      Node2 Element2 g2 H2 elind2 Eind2 zs2 NN2 sN2...
                      g3 H3 zs3 NN3 sN3  g4 H4 zs4 NN4 sN4






