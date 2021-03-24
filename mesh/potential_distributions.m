% Plotting potential distributions


load Ucur.mat


ng = max(HF1st(:));
U_dis = Ucur(1:ng,3); 


figure(1)

%%%%%%%%%%%%%%%%%%%%%%%
%%% Vertical slices %%%
%%%%%%%%%%%%%%%%%%%%%%%

IntpMatrices2 = []; % set IntpMatrices empty only in the beginning. 
                    % Drawslices.m returns the IntpMatrices; 
                    % after this drawing the slices happens quickly
                    
xyslices2 = max(gF1st(:,3))*[0:.4:1]; 
xzslices2 = max(gF1st(:,2))*[0];
yzslices2 = max(gF1st(:,1))*[.5];


viewangle = [0 0];
IntpMatrices2 = DrawSlices(gF1st,HF1st,NodeF1st,U_dis,xyslices2,xzslices2,yzslices2,'cylinder',max(gF1st(:,1))/50,101,IntpMatrices2,viewangle);


simp_plot_3dvec(gF1st,HF1st,'p(:,3)<5',U_dis,U_dis);

%%%%%%%%


tht = -pi/4;

Rot_z = [cos(tht) -sin(tht) 0; sin(tht) cos(tht) 0; 0 0 1]; 
Nnodes=max(size(gF1st));

gtemp = zeros(size(gF1st));



for ii = 1:Nnodes
      
    gtemp(ii,:) =(Rot_z*gF1st(ii,:)')'; 
    
end

figure(1),simp_plot_3d(gF1st,HF1st), view(40,12)
figure(2),simp_plot_3d(gtemp,HF1st), view(40,12)


IntpMatrices2 = []; % set IntpMatrices empty only in the beginning. 
                    % Drawslices.m returns the IntpMatrices; 
                    % after this drawing the slices happens quickly
                    
xyslices2 = max(gtemp(:,3))*[0:.4:1]; 
xzslices2 = max(gtemp(:,2))*[0];
yzslices2 = max(gtemp(:,1))*[.5];


viewangle = [0 0];
%IntpMatrices2 = DrawSlices(ginv,Hinv,Nodeinv_P,sigma,xyslices2,xzslices2,yzslices2,'cylinder',max(ginv(:,1))/50,101,IntpMatrices2,viewangle);
IntpMatrices2 = DrawSlices(gtemp,HF1st,NodeF1st,U_dis,xyslices2,xzslices2,yzslices2,'cylinder',max(gtemp(:,1))/50,101,IntpMatrices2,viewangle);





% %simp_plot_3d(g_chopped,H_chopped)
% simp_plot_3dvec(gF1st,HF1st,'p(:,3)<20',U_dis,U_dis);
% 
% IntpMatricesF = [];
% xyslicesF = max(gF1st(:,3))*[0:.1:1]; 
% xzslicesF = max(gF1st(:,2))*[];
% yzslicesF = max(gF1st(:,1))*[];
% 
% 
% 
% [VisElementsF, VisNodesF, IntpMatricesF] = Generate2DSlices(gF1st,HF1st,NodeF1st,U_dis,xyslicesF,xzslicesF,yzslicesF,max(gF1st(:,1))/50,IntpMatricesF);
% 
% PlotMyResults(VisElementsF, VisNodesF, IntpMatricesF, U_dis,[min(U_dis), max(U_dis)],101)
% 
% 
% 
% IntpMatrices2 = []; % set IntpMatrices empty only in the beginning. 
%                     % Drawslices.m returns the IntpMatrices; 
%                     % after this drawing the slices happens quickly
% 
% 
% 
% viewangle = [35 40];
% [VisElementsF, VisNodesF, IntpMatricesF] = Generate2DSlices(gF1st,HF1st,NodeF1st,U_dis,xyslicesF,xzslicesF,yzslicesF,max(gF1st(:,1))/50,IntpMatricesF);



% IntpMatricesF = [];
% xyslicesF = max(g1st(:,3))*[0:.1:1]; 
% xzslicesF = max(g1st(:,2))*[];
% yzslicesF = max(g1st(:,1))*[];
% 
% range= [min_sigma, max_sigma];
% 
% [VisElementsF, VisNodesF, IntpMatricesF] = Generate2DSlices(gF1st,HF1st,NodeF1st,theta_true,xyslicesF,xzslicesF,yzslicesF,max(gF1st(:,1))/50,IntpMatricesF);
% PlotMyResults(VisElementsF, VisNodesF, IntpMatricesF, theta_true,range,101)
% 
% 
% IntpMatrices = [];
% xyslices = max(ginv(:,3))*[0:.1:1]; 
% xzslices = max(ginv(:,2))*[];
% yzslices = max(ginv(:,1))*[];
% 
% [VisElements, VisNodes, IntpMatrices] = Generate2DSlices(ginv,Hinv,Nodeinv_P,theta(:,end),xyslices,xzslices,yzslices,max(ginv(:,1))/50,IntpMatrices);
% 
% 
% PlotMyResults(VisElements, VisNodes, IntpMatrices, theta(:,end),range,101)
% 
