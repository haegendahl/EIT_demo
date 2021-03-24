function [H,g,h,Node,Element,elind,eltetra,E] = makecylmeshHalf(GEOMETRYTYPE,geofile,order,ng_parameters,meshname,MaxElemSizeUnderElectrode)


%keyboard

if nargin<3
  order = 1;
end

if ~((order==1) || (order==2)) error('Unsupported mesh order!'), end

% run geo-file construction script
if strcmp(GEOMETRYTYPE,'kiekko_circ_electodes');
  Nel = 16;
  [elc,elp] = CreateCircMesh(geofile,[1:Nel]);
elseif strcmp(GEOMETRYTYPE,'tankki_square_electodes2rows');
  % TAMA PITAA MUUTTAA:
  CreateCircMesh_2DbarmultibarsTankki2rows(geofile,ms);
elseif strcmp(GEOMETRYTYPE,'watertank');
  Nel = 16;
  %CreateFlatTankMeshNobars(geofile);
  CreateFlatTankMeshNobarsDenseElectrodes(geofile,MaxElemSizeUnderElectrode);
end


% define geometry and mesh files for input & output
%geofile = 'planemesh_1bar_test2.geo';
%meshfile = 'barmesh.vol';

% set meshing options for Netgen (in a cell-array)
% These are the only supported parameters so far...
% (Note: names are case insensitive but extra spaces are not allowed)


% run Netgen
start_netgen(geofile,meshname,ng_parameters);

% read mesh and build structures
[H,g,h,tmp] = read_ng(meshname,order);
H = H';
g = g';
h = h';
surf = h(:,1);
bcnr = h(:,2);
if(order==1)
  h = h(:,3:5);
elseif(order==2)
  h = h(:,3:8);
end

%figure,simp_plot_3d(g,H)


DRAWELECTRODES = 1;

if DRAWELECTRODES
  figure(1),clf
  hold on
  for ii = 100+(1:max(bcnr)-1)
   I2 = find(bcnr==ii);
   patch('faces',h(I2,:),'vertices',g,'facecolor',[1 0 0]),view(3), axis equal
   disp(ii)
   pause
  end
end


% get elind
elind = cell(Nel/2,1); 
if strcmp(GEOMETRYTYPE,'kiekko_circ_electodes');
  for ii=1:Nel/2
    I = find(bcnr==ii);
    f = h(I,:);
    elind{ii} = unique(f(:));
  end
elseif strcmp(GEOMETRYTYPE,'watertank');
  for ii=1:Nel/2
    I = find(bcnr==100+ii);
    f = h(I,:);
    elind{ii} = unique(f(:));
  end
end

Node = MakeNode3dSmallFast(H,g);
[eltetra,E] = FindElectrodeElements2(H,Node,elind,order);
%[eltetra,E] = FindElectrodeElements(H,elind,Nel,order);
Element = MakeElement3dSmallCellFast(H,eltetra,E);

eval(['save ',meshname,' g H Node Element'])

%figure(100),clf,hold on
%for iii = 1:Nel
%  iii
%  ii = elind{iii};
%  plot3(g(ii,1),g(ii,2),g(ii,3),'+')
%  axis image
%  pause
%end


