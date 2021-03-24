function DrawTarget(ginv,Hinv,Nodeinv_P,theta,fig1,fig2, viewangle,cmin,cmax)

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Horizontal slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

IntpMatrices = []; % set IntpMatrices empty only in the beginning. 
                   % Drawslices.m returns the IntpMatrices; 
                   % after this drawing the slices happens quickly

xyslices = max(ginv(:,3))*[0]; 
xzslices = max(ginv(:,2))*[0];
yzslices = max(ginv(:,1))*[0];

%viewangle = [35 40];
IntpMatrices = DrawSlices(ginv,Hinv,Nodeinv_P,theta(:,end),xyslices,xzslices,yzslices,'cylinder',max(ginv(:,1))/50,fig1,IntpMatrices,viewangle);
caxis([cmin,cmax])


%%%%%%%%%%%%%%%%%%%%%%%
%%% Vertical slices %%%
%%%%%%%%%%%%%%%%%%%%%%%

IntpMatrices2 = []; % set IntpMatrices empty only in the beginning. 
                    % Drawslices.m returns the IntpMatrices; 
                    % after this drawing the slices happens quickly

xyslices2 = max(ginv(:,3))*[0:.1:1]; 
xzslices2 = max(ginv(:,2))*[];
yzslices2 = max(ginv(:,1))*[];

%viewangle = [35 40];
IntpMatrices2 = DrawSlices(ginv,Hinv,Nodeinv_P,theta(:,end),xyslices2,xzslices2,yzslices2,'cylinder',max(ginv(:,1))/50,fig2,IntpMatrices2,viewangle);
caxis([cmin,cmax])


%


