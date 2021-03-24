function [ginv,Hinv] = MakeUniformInvGrid(cylinder_height,cylinder_radius,radius_div,drawmode,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh for inverse problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r = 7.5*[1 .87 .69 .52 .36  0.18,  0];
%NN = [48 27 24 21 14 6 1];

%r = 7.5*[1 .88 .79 .66 .52 .36 .22 .14,  0];
%NN = [64  48 44 36 24 20 16 10 1];

r = cylinder_radius/radius_div*[radius_div:-1:0];
elem_width = abs(diff(r(1:2)));
NN = round(2*pi*r/elem_width);
NN(1) = NN(1)-mod(NN(1),2);
NN(end) = 1;

Nlayers = round(cylinder_height/elem_width);
zs = cylinder_height/Nlayers*ones(Nlayers,1);

el = [1,1];
L = NN(1)/2;
 
h = sum(zs);
N = zeros(1,Nlayers);
N(round(Nlayers/2)) = 1;
%el_width = el(1)*2*pi*max(r)/(L*sum(el));
%bel_width = el(2)*2*pi*max(r)/(L*sum(el));

%[Eindtmp,elindtmp,ginv,Hinv] = MakeMesh3D_FlowLoop(N,r,NN,h,el,zs,L,el_width,bel_width);

[Eindtmp,elindtmp,ginv,Hinv] = MakeMesh3D(N,r,NN,h,el,zs);


if(strcmp(drawmode,'draw'))
 figure,
 trimesh(Hinv,ginv(:,1),ginv(:,2),ones(size(ginv,1),1),ones(size(ginv,1),1)),axis image
 view(2)
 %DrawCylindricalGrid(ginv,Hinv,elindtmp)
 figure(101),simp_plot_3d(ginv,Hinv)
end


eval(['save ',filename,' ginv Hinv'])



