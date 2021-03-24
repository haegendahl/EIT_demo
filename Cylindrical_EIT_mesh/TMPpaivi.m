path(path,'../../../EIT_codes/mesh')


geofile = 'tmpgeo.geo';
meshname = 'tmpmesh'
radius = 10;
height = 5;
order = 1;

ng_parameters = {'Mesh granularity', 2,...   % overall mesh density, 1 = coarsest, 5 = finest
                 '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
                 'Max mesh-size', 5,...   % maximum allowed element size
                 'Mesh-size grading', 0.2,... 
                 'Elements per curvature', 1.5,... 
                 'Elements per edge', 2.0};
maxH = 0.5;



MakeGeofile_circ_tmp(geofile,radius,height,maxH)
% run Netgen
start_netgen(geofile,meshname,ng_parameters);
[H,g,h,tmp] = read_ng(meshname,order);

H = H';
g = g';

simp_plot_3d(g,H)


