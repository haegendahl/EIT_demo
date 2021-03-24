 function M=ReadVolMesh2ndOrder(varargin)
 
 if nargin ==0,
     [filename, pathname] = uigetfile('*.vol', 'Select mesh file','mesh');
 else
     filename=varargin{1};
     pathname=[];
 end
 %keyboard
 drawnow
 fid=fopen([pathname,filename]);
 
 LINE='';
 while ~strcmp(LINE,'surfaceelementsgi'),
     LINE=fgetl(fid);
 end
 
 Nsurfelem=str2num(fgetl(fid));
 SURFMAT=zeros(Nsurfelem,17);
 for ii=1:Nsurfelem,
     SURFMAT(ii,:)=str2num(fgetl(fid));
 end

 LINE='';
 while ~strcmp(LINE,'volumeelements'),
     LINE=fgetl(fid);
 end
 Nelem=str2num(fgetl(fid));
 VOLMAT=zeros(Nelem,12);
 for ii=1:Nelem,
     VOLMAT(ii,:)=str2num(fgetl(fid));
 end

 LINE='';
 while ~strcmp(LINE,'points'),
     LINE=fgetl(fid);
 end
 Nnode=str2num(fgetl(fid));
 G=zeros(Nnode,3);
 for ii=1:Nnode,
     G(ii,:)=str2num(fgetl(fid));
 end
 
% LINE='';
% while ~strcmp(LINE,'identifications'),
%   LINE=fgetl(fid);
% end
% Nidentifications=str2num(fgetl(fid));

% identifications=zeros(Nidentifications,3);
% for ii=1:Nidentifications,
%     identifications(ii,:)=str2num(fgetl(fid));
% end
 
 fclose(fid);
 M.g=G;
 M.t=VOLMAT;
 M.e=SURFMAT;
 M.Nnode=Nnode;
 M.Nelem=Nelem;
 M.Nsurf_elem=Nsurfelem;
 M.fname=filename;
% M.identifications = identifications;
