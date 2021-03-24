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

zs = el_width_z/2*ones(2*(2*4+4+3*2),1);
r = L*(el_width+bel_width)/(2*pi)*[1 .87 .69 .52 .36  0.18,  0];
NN = [48 27 24 21 14 6 1];

el = [2,1]; 
h = sum(zs);
N = [0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0  0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0];

[Eind,elind,g,H] = MakeMesh3D_FlowLoop(N,r,NN,h,el,zs,L,el_width,bel_width);

%figure,
%trimesh(H,g(:,1),g(:,2),ones(size(g,1),1),ones(size(g,1),1)),axis image
%view(2)

sN = max(size(g));

Draw_Cylindrical_Grid

[Element]=MakeElement3dSmall(H,elind,Eind);
[Node]=MakeNode3dSmall(Element,g);

sN = max(size(g));
WW = speye(sN,sN);
Agrad = jacobian3dnode(Node,Element,WW);


%eval(['save ',meshfile,' Node Element g H N elind Eind Agrad zs NN sN'])       



%%% Grid for 2nd order basis %%%


if (2nd_order_basis == 'yes')

  [g3,H3,elind3,Eind3] = mesh2nd(g,H,elind,Eind,L);

  elindb3=elind3;
  Eindb3=Eind3;
  for jj=1:L   
   eval(['elind3.ElecNo', int2str(jj),'=elindb3(jj,:);']) 
   eval(['Eind3.ElecNo', int2str(jj),'=Eindb3(jj,:);'])
  end

  eval(['save ',meshfile,' Node Element g H N elind Eind Agrad zs NN sN g3 H3 elind3 Eind3'])


else

   eval(['save ',meshfile,' Node Element g H N elind Eind Agrad zs NN sN'])

end

