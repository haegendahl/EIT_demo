function [] = CreateFlatTankMeshNobarsDenseElectrodes(geofile,MaxElemSizeBoundEl)
% 
% sylinterihila, suorakaide-elektrodit
%

%keyboard

elec_s = 2.5;  % electrodin leveys
elec_h = 7;  % electrodin korkeus
R = 28/2;    % hilan säde
elec_v = 4;  % electroditasojen väli (yläreunasta alareunaan)
taso_l = 1;  % electroditasojen lukumäärä
elec_l = 16; % electrodien lukum. tasolla
tila_a = 0;  % tyhjä tila hilan pohjalta elektrodiin
tila_y = 0;  % tyhjä tila hilan pinnalta elektrodiin
radius = R;

% concrete disc data (cm)
%radius = 7.5;  
%height = 3;
%Nel = 16;
%elr = 0.5;  % electrode radius
%z = height/2;
%z0 = num2str(z);
maxh = 1;

% calculate boundary electrode centers
% (and corresponding cylinder end points for netgen)
%tht = linspace(0,2*pi,Nel+1); tht(Nel+1) = [];
%rxy = diag([1.2*radius 0.8*radius]);
%cth = cos(tht(:)); sth = sin(tht(:));
%elcenters = radius*[cth sth];
%elpoints = [[cth cth]*rxy [sth sth]*rxy];

Nel=taso_l*elec_l; % electrodien lkm.
bcnr = 1:Nel;
elmaxh=MaxElemSizeBoundEl*ones(Nel,1);

l=5; % laatikon syvyys
Thta=zeros(elec_l,1);
Thta(1)=0;
for ii=2:elec_l,
    Thta(ii)=Thta(ii-1)+2*pi/(elec_l);
end

height = tila_a + tila_y + taso_l*elec_h + (taso_l-1)*elec_v;


elec_th=elec_s/(R); % electrodin leveys radiaaneina

xx=zeros(2,elec_l);
yy=zeros(2,elec_l);

for ii=1:elec_l,
    thta1=Thta(ii)-elec_th/2;
    thta2=Thta(ii)+elec_th/2;
    thta=[thta1 thta2];
    [x y]=pol2cart(thta,R*ones(1,2));
    xx(:,ii)=x(:);
    yy(:,ii)=y(:); 
end

xx=repmat(xx,1,taso_l);
yy=repmat(yy,1,taso_l);

z1=tila_a+elec_h;
z2=z1*ones(2,elec_l);
zz=[];
zz=[zz z2];
z2=[];

for ii=2:taso_l,
    z1=z1+elec_v+elec_h;
    z2=z1*ones(2,elec_l);
    zz=[zz z2];
end

XX=xx;
YY=yy;
ZZ=zz;


thta_rep=repmat(Thta,1,taso_l);

pisteet=[];
%alpha = 0.2105;

%figure, hold on


for ii=1:Nel,
    
     
     xx=XX(:,ii);
     yy=YY(:,ii);
     z=ZZ(:,ii);
     z(2)=z(1)-elec_h;


    Ett=[];
    for jj=1:max(size(xx)),
        ett=sqrt((xx-xx(jj)).^2+(yy-yy(jj)).^2);
        Ett=[Ett ett];
    end
    
    d=max(max(Ett))
    
    %plot3(x,y,z,'*')
    h=max(z)-min(z);

    alpha=elec_th;
    
    %alpha=max(thta1)-min(thta1);
    X=tan(0.5*alpha)*(R-0.5*l);
    Y=sqrt(X^2+(R-0.5*l)^2);
    Z=0.5*d-X;
    r1=sqrt((R-0.5*l)^2+(0.5*d)^2);
    beta=acos(-(Z^2-r1^2-Y^2)/(2*r1*Y));
    
    r2=sqrt((R+0.5*l)^2+(0.5*d)^2);
    beta2=asin((0.5*d)/(R+0.5*l));
    

    tht1=thta_rep(ii)+0.5*alpha+beta;
    
    if (tht1>2*pi)
        tht1=tht1-2*pi;
    end
    tht5=tht1;
    
    %tht4=thta-0.5*alpha-beta;
    tht4=thta_rep(ii)-0.5*alpha-beta;
    if (tht4<0)
        tht4=2*pi-abs(tht4);
    end
    tht8=tht4;
   
    %tht2=thta+beta2;
    tht2=thta_rep(ii)+beta2;
    if (tht2>2*pi)
        tht2=tht2-2*pi;
    end
    tht6=tht2;
    
    %tht3=thta-beta2;
    tht3=thta_rep(ii)-beta2;
    if (tht3<0)
        tht3=2*pi-abs(tht3);
    end
    tht7=tht3;
    
 
    
    r4=r1;
    r3=r2;
    r5=r1;
    r6=r2;
    r7=r3;
    r8=r4;
    
    z1=max(z);
    z2=z1;
    z3=z1;
    z4=z1;
    
    z5=min(z);
    z6=z5;
    z7=z5;
    z8=z5;
    
    % valitaan 1 ja 7
    
    if tht1>pi
        tht1=-pi+(tht1-pi);
    end
    
    if tht7>pi
        tht7=-pi+(tht7-pi);
    end
    
    [xx1 yy1 zz1]=pol2cart(tht1,r1,z1);
    [xx2 yy2 zz2]=pol2cart(tht2,r2,z2);
    [xx3 yy3 zz3]=pol2cart(tht3,r3,z3);
    [xx4 yy4 zz4]=pol2cart(tht4,r4,z4);
    [xx5 yy5 zz5]=pol2cart(tht5,r5,z5);
    [xx6 yy6 zz6]=pol2cart(tht6,r6,z6);
    [xx7 yy7 zz7]=pol2cart(tht7,r7,z7);
    [xx8 yy8 zz8]=pol2cart(tht8,r8,z8);
    
  
    
    xxx=[xx1 xx2 xx3 xx4 xx5 xx6 xx7 xx8];
    yyy=[yy1 yy2 yy3 yy4 yy5 yy6 yy7 yy8];
    zzz=[zz1 zz2 zz3 zz4 zz5 zz6 zz7 zz8];
    %plot3(xxx,yyy,zzz,'-*')
    
    pist=[xxx(:) yyy(:) zzz(:)];
    pist=pist';
    pisteet=[pisteet pist];

end

cameramenu     
%axis equal

[boxstr] = asavec2str_B(pisteet',Nel);

elstr = [];
for ii=1:Nel
    elstr = [elstr, sprintf('solid el%d = box%d and cyl;\n',ii,ii)];
end

tlostr = [];
for ii=1:Nel
    tlostr = [tlostr, ...
              sprintf('tlo el%d cyl -col=[1,0,0] -bc=%d;\n',ii,ii);];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nintel stuff

%elmaxh=.5*ones(Nel,1);
%barmaxh = 0.25;
%% steel bar data
%bx = meshstruct.barxy(:,1);
%by = meshstruct.barxy(:,2);
%br = meshstruct.barr;
%Nintel = meshstruct.Nintel;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write into geo-file... %
%%%%%%%%%%%%%%%%%%%%%%%%%%

elstr = [];
for ii=1:Nel
    elstr = [elstr, sprintf('solid el%d = box%d and cylwall;\n',ii,ii)];
end

tlostr = [];
for ii=1:Nel
    tlostr = [tlostr, ...
     sprintf(['tlo el%d cylwall -col=[1,0,0] -bc=%d -maxh=',num2str(elmaxh(ii)),';\n'],ii,ii+100);];
end


% open file
fid = fopen(geofile,'wt');  % open for writing (text mode)
if fid==-1
  error('Cannot open file: %s',fname)
else
  fprintf(fid,'algebraic3d\n\n#box definitions\n');
  
  for ii=1:Nel
    fprintf(fid,'%s\n',boxstr{ii});
  end
  
  fprintf(fid,['\nsolid cylwall = cylinder(0,0,0; 0,0,',num2str(height),'; %0.4g);\n'],R);
  fprintf(fid,['solid main = cylwall and plane(0,0,0;0,0,-1) and plane(0,0,',num2str(height),';0,0,1);\n']);
  
  fprintf(fid,'\n%s\n',elstr);
  
  fprintf(fid,'\n# top level objects\n');
  fprintf(fid,'tlo main -transparent;\n\n');
  fprintf(fid,'%s\n',tlostr);
  fclose(fid);
end


