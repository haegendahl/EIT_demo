
close all
clear all
MakePaths

%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%

GEOMETRYTYPE = 'watertank3D_2'
Nel = 64;
   
FORWARD_MESH_WITH_TARGET_READY = 0;
fwd1stmeshname = 'TANK_MESH_1st';
fwd2ndmeshname = 'TANK_MESH_2nd';
%geofilename = 'kiekkohila.geo';
geofilename = 'watertank3D_2.geo';
%agradfilename = 'kiekkohilaharvaAgrad';
agradfilename = 'watertank3DAgrad';
order = 2;

ng_parameters = {'Mesh granularity', 5,...   % overall mesh density, 1 = coarsest, 5 = finest
                 '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
                 'Max mesh-size', 2,...   % maximum allowed element size
                 'Mesh-size grading', .2,... 
                 'Elements per curvature', 2,... 
                 'Elements per edge', 2.0};
MaxElemSizeUnderElectrode = .5;
%MaxElemSizeUnderElectrode = .25;


INVMESHREADY = 1;
INTMATRIX_INV_1ST_READY = 1;
%invmeshtype = 'Forwardmesh';
invmeshtype = 'Uniform';
elementsize_inv = 1.5;
invmeshfile = 'TANK_UNIFORM_INV';
intpmatfile_inv_1st = 'Mapping_to_unimesh_1st';

%%%%%%%%%%%%%%%%%%
% Forward meshes %
%%%%%%%%%%%%%%%%%%

if ~FORWARD_MESH_WITH_TARGET_READY

  [H,g,h2nd,Node,Element,elind2nd,eltetra2nd,E2nd] = makecylmesh(GEOMETRYTYPE,geofilename,2,ng_parameters,fwd2ndmeshname,MaxElemSizeUnderElectrode,Nel);
  [H1st,g1st,Node1st,Element1st,elind1st,eltetra1st,E1st] = Reduce2ndOrderMesh(H,g,elind2nd,fwd1stmeshname);
  figure(2),simp_plot_3d(g1st,H1st), view(40,12)
  
  %% For the Jacobian
  WW = speye(size(g1st,1));
  Agrad = jacobian3dnodeFast(Node1st,Element1st,WW);
  eval(['save ',agradfilename,' Agrad'])

end

figure(3),clf,simp_plot_3d(g,H)