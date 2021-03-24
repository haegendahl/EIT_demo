function [volmesh,fem]=GetMesh_netgen2ndOrder(varargin)


  if nargin == 0
    M = ReadVolMesh2ndOrder;
  else
    M = ReadVolMesh2ndOrder(varargin{1});
  end
  
  
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% First mesh dependent stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % mesh consist of one volume, next extract indices to certain domain
  
  dms = M.t(:,1);
  i1 = find(dms==1);
 
  volmesh.domains=zeros(size(M.t(:,1)));
  volmesh.domains(i1)=1;
   
  volmesh.g=M.g;
  volmesh.H=M.t(:,3:end);
    
  %% surfaces
  srf = M.e(:,1);
  srf_numbers = unique(srf);
  
  volmesh.surfaces=zeros(size(M.e(:,1)));
  
  for kk = 1:size(srf_numbers,1)
    s1 = find(srf==srf_numbers(kk));
    volmesh.surfaces(s1)=kk;
  end
  
  volmesh.E=M.e(:,6:11);
  volmesh.Nlayer=max(volmesh.domains); 
  
  %for ii=1:max(volmesh.surfaces),
  %  xyz1=volmesh.g(volmesh.E(find(volmesh.surfaces==ii),:),:);
  %  [a,b,rr]=cart2sph(xyz1(:,1),xyz1(:,2),xyz1(:,3));
  %  RR(ii)=mean(rr);
  %end

  % some dimensions
  volmesh.Nelem=M.Nelem;
  volmesh.Nnode=M.Nnode;
  volmesh.fname=M.fname;
  
  fem.mesh.t=[volmesh.H';volmesh.domains'];
  fem.mesh.e=[volmesh.E zeros(size(volmesh.E,1),9)]';
  fem.mesh.e(10,:)=volmesh.surfaces;
  fem.mesh.p=volmesh.g';
  fem.mesh.h=M.e;
  fem.geom=[];
 
% $$$  figure(1), clf,
% $$$  meshplot(fem,'ellogic','y>0'), axis equal;
% $$$  set(gca,'zdir','reverse')
% $$$  
% $$$  %figure(2),clf,
% $$$  %meshplot(fem,'ellogic','y>2','bdl',[1]), axis equal;
% $$$  %set(gca,'zdir','reverse')
% $$$  
% $$$  %figure(3),clf,
% $$$  %meshplot(fem,'ellogic','y>2','BoundMode','on'), axis equal;
% $$$  %set(gca,'zdir','reverse')
% $$$  
% $$$  figure(4),clf,
% $$$  meshplot(fem,'edgemode','off','boundmode','off','dboundmode','on','curvemode','on'), axis equal;
% $$$  set(gca,'zdir','reverse')
