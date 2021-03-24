%%% Mesh for Forward simulations %%%

zs2 = .5*ones(36,1);
r2 = 5.5*[1, 0.85, 0.65, 0.45, 0.25, 0];
NN2 = [64 48 24 12 6 1];
el2 = [2,2]; 
h2 = sum(zs2);
N2 = [0 0 0 0  0 0 0 0 1  1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0];

[Eind2,elind2,g2,H2] = MakeMesh3D(N2,r2,NN2,h2,el2,zs2);

figure,
trimesh(H2,g2(:,1),g2(:,2),ones(size(g2,1),1),ones(size(g2,1),1)),axis image
view(2)

%[Element2]=MakeElement3dSmall(H2,elind2,Eind2);
%[Node2]=MakeNode3dSmall(Element2,g2);

sN2 = max(size(g2));
