function [] = MakeGeofile_circ_tmp(geofile,radius,height,maxH)
% [elcenters] = CreateCircularMesh(geofile,meshstruct)
%    meshstruct = struct('radius',
 

% concrete disc data (cm)
%radius = 7.5;  
%height = 3;
%Nel = 16;
%elr = 0.5;  % electrode radius
z = height/2;
z0 = num2str(z);
maxh = 1;
 
% calculate electrode centers
% (and corresponding cylinder end points for netgen)
%tht = linspace(0,2*pi,Nel+1); tht(Nel+1) = [];
%rxy = diag([1.2*radius 0.8*radius]);
%cth = cos(tht(:)); sth = sin(tht(:));
%elcenters = radius*[cth sth];
%elpoints = [[cth cth]*rxy [sth sth]*rxy];

%elmaxh=.5*ones(Nel,1);
%barmaxh = 0.25;


% open geofile and write the data
ngfid = fopen(geofile,'w');
if ngfid==-1, error('couldn''t open geofile!'), end


% header
fprintf(ngfid,'algebraic3d\n\n');

Nel = 0;

% main object
fprintf(ngfid,'solid cyl = cylinder(0,0,0; 0,0,%0.3f; %0.3f) -maxh=%0.3f;\n',height,radius,maxH);
fprintf(ngfid,'solid bigcyl = plane(0,0,0; 0,0,-1)\n\t and ');
fprintf(ngfid,'plane(0,0,%0.3f; 0,0,1)\n\t and cyl -bc=%d;\n\n',...
        height,Nel+2);


% electrodes
%for ii=1:Nel,
%  x1 = elpoints(ii,1); y1 = elpoints(ii,3);
%  x2 = elpoints(ii,2); y2 = elpoints(ii,4);
%  n1 = [x1 y1]-[x2 y2];  n1 = n1/norm(n1);
%  n2 = -n1;
%  pl1str = [num2str(x1),',',num2str(y1),',',z0,';',num2str(n1(1)),','...
%            num2str(n1(2)),',0'];
%  pl2str = [num2str(x2),',',num2str(y2),',',z0,';',num2str(n2(1)),','...
%            num2str(n2(2)),',0'];  
%  cylstr = [num2str(x2),',',num2str(y2),',',z0,';',num2str(x1),','...
%            num2str(y1),',',z0,';',num2str(elr)];
%  
%  fprintf(ngfid,'solid cyl%d = plane(%s)\n\t and plane(%s)\n\t and cylinder(%s);\n',...
%          ii,pl1str,pl2str,cylstr);
%end
%
%for ii=1:Nel,
%  fprintf(ngfid,'solid el%d = bigcyl and cyl%d;\n',ii,ii);
%end
%
%% top level objects
%fprintf(ngfid,'\ntlo bigcyl -maxh=%0.1f;\n',maxh);
%for ii=1:Nel,
%  fprintf(ngfid,'tlo el%d cyl -maxh=%0.1f -col=[1,0,0] -bc=%d;\n',...
%          ii,elmaxh(ii),bcnr(ii));
%end
%fprintf(ngfid,'tlo bar barel -bc=%d;\n',Nel+1);

fprintf(ngfid,'tlo bigcyl; \n');

fclose(ngfid);

