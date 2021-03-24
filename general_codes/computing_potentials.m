
        
        %------- Computing the inner potentials
        
        clc

        clear all
        
        close all
        
        MakePaths
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Mesh parameters %%%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % - FORWARD MESH for Simulating the EIT Measurements
        
        GEOMETRYTYPE = 'watertank3D_2'
        %geofilename = 'watertank3D_2.geo'; % sphere inclusion
        geofilename = 'watertank3D_3.geo'; % cylinder inclusion
        %geofilename = 'wt_triangular_prism.geo'; % triangular prism
       
        %geofilename = 'watertank3D.geo';
        
        Nel = 64;
        
        FORWARDMESHREADY = 0;
        fwd1stmeshname = 'TANK_MESH_1st';
        fwd2ndmeshname = 'TANK_MESH_2nd';
        
        order = 2;
        
        ng_parameters = {'Mesh granularity', 5,...   % overall mesh density, 1 = coarsest, 5 = finest
            '2nd order elements', order-1,... % 0 or 1 (1 = use 2nd order)
            'Max mesh-size', 1,...   % maximum allowed element size
            'Mesh-size grading', .2,...
            'Elements per curvature', 2,...
            'Elements per edge', 2.0};
        MaxElemSizeUnderElectrode = .5;
        
        
        %%%%%%%%%%%%%%%%%%
        % Forward meshes %
        %%%%%%%%%%%%%%%%%%
        
        if ~FORWARDMESHREADY
            
            [HF,gF,h2ndF,NodeF,ElementF,elind2ndF,eltetra2ndF,E2ndF] = makecylmesh(GEOMETRYTYPE, geofilename,2,ng_parameters,fwd2ndmeshname,MaxElemSizeUnderElectrode,Nel);
            [HF1st,gF1st,NodeF1st,ElementF1st,elindF1st,eltetraF1st,EF1st] = Reduce2ndOrderMesh(HF,gF,elind2ndF,fwd1stmeshname);
          %  figure(100),simp_plot_3d(gF1st,HF1st), view(40,12)
            
            %% For the Jacobian
            %   WW = speye(size(gF1st,1));
            %   Agrad = jacobian3dnodeFast(Node1st,Element1st,WW);
            %   eval(['save ',agradfilename,' Agrad'])
            
        else
            
            mS = load(fwd1stmeshname);
            gF1st = mS.g;
            HF1st = mS.H;
            NodeF1st = mS.Node;
            ElementF1st = mS.Element;
            eltetraF1st = mS.eltetra;
            EF1st = mS.E;
            clear mS
            
            
            mS = load(fwd2ndmeshname);
            gF = mS.g;
            HF = mS.H;
            NodeF = mS.Node;
            ElementF = mS.Element;
            clear mS
            
        end
        
        % 2nd order organizing...
        [ElementF] = organize_faces2nd(ElementF,gF);
        [HF ElementF] = organize_tetras2nd(HF,ElementF,gF);
       
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GENERATING THE TARGET %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        min_sigma = 2;
        max_sigma = 2;
        
        
        %cylinder
        x0=5;
        y0=3;
        r=3;
        index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2 <= (r+1e-2)^2;
       
        theta_true = max_sigma*ones(size(gF1st,1),1);
        theta_true(index) = min_sigma;
        
        
        %----
        
        figure(1)
        set(gcf,'renderer','zbuffer')
        [hp1,hp2] = simp_plot_3dvec(gF1st,HF1st,'p(:,3)<5',theta_true,theta_true);
        view(40,12)
        axis tight
        colorbar
     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Simulate EIT measurements %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        L = 64                                        % Number of electrodes
        % Number of electrodes
        z = 0.00001*ones(L,1);          % Contact Impedances
        
        %disp(['Using: Opposite current injections and measurements',10]);
        % I = toeplitz([1; zeros(L/2-1,1); -1; zeros(L/2-1,1)]',[1 zeros(1,L/2-1)]); % Opposite current injections
        disp(['Using: Custom current injections and adj measurement pattern',10]);
        I1 = -eye(L/4,L/4);
        
        I2 = zeros(L/2,L/4);
        
        I3 = toeplitz([zeros(1,(L/8)), 1  zeros(1,L/8-1)]);
        
        I = [I1;I2;I3]; % First layer oppsite last layer
        
        %MeasPatt = [ones(1,L-1);-eye(L-1)];                                        % Opposite measurements
        
        MeasPatt = Current(64,0,'adj');
        
        U = ForwardSolution3d2ndElectrode_G(NodeF,ElementF,I,theta_true,z,MeasPatt,'real');
        Uel= U.Electrode;
        Umes = U.MeasField;
        Ucur = U.Current;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Ploting Potential Distribution  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(103)
        Uel = Uel(:);
        plot(Uel,'r')
        
        ng = max(HF1st(:));
        U_dis = Ucur(1:ng,7);
  
        
        tht = -3*pi/4;
        
        Rot_z = [cos(tht) -sin(tht) 0; sin(tht) cos(tht) 0; 0 0 1];
        Nnodes=max(size(gF1st));
        
        gtemp = zeros(size(gF1st));
        
                
        for ii = 1:Nnodes
            
            gtemp(ii,:) =(Rot_z*gF1st(ii,:)')';
            
        end
        
        figure(1),simp_plot_3d(gF1st,HF1st), view(40,12)
        figure(2),simp_plot_3d(gtemp,HF1st), view(40,12)
        
       % simp_plot_3dvec(gtemp,HF1st,'p(:,3)<5',theta_true,theta_true)
        %---- Generating Slices
        IntpMatrices2 = []; % set IntpMatrices empty only in the beginning.
        % Drawslices.m returns the IntpMatrices;
        % after this drawing the slices happens quickly
        
        xyslices2 = max(gtemp(:,3))*[0:.4:1];
        xzslices2 = max(gtemp(:,2))*[0];
        yzslices2 = max(gtemp(:,1))*[.5];
        
        figure(101)
        viewangle = [0 0];
        IntpMatrices2 = DrawSlices(gtemp,HF1st,NodeF1st,U_dis,xyslices2,xzslices2,yzslices2,'cylinder',max(gtemp(:,1))/50,2,IntpMatrices2,viewangle);
        colorbar
        title('without inclusion')
       
  %--------------- Second part      
     
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GENERATING THE TARGET %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        min_sigma = 0.5;
        max_sigma = 2;
        
        
        %cylinder
        x0=-5;
        y0=5;
        r=4;
        index=(gF1st(:,1)-x0).^2+(gF1st(:,2)-y0).^2 <= (r+1e-2)^2;
       
        theta_true = max_sigma*ones(size(gF1st,1),1);
        theta_true(index) = min_sigma;
        
        
        %----
        
        figure(101)
        set(gcf,'renderer','zbuffer')
        [hp1,hp2] = simp_plot_3dvec(gF1st,HF1st,'p(:,3)<21',theta_true,theta_true);
        view(40,12)
        axis tight
        colorbar
     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Simulate EIT measurements for second part %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        U = ForwardSolution3d2ndElectrode_G(NodeF,ElementF,I,theta_true,z,MeasPatt,'real');
        Uel= U.Electrode;
        Umes = U.MeasField;
        Ucur = U.Current;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Ploting Potential Distribution  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(104)
        Uel = Uel(:);
        plot(Uel,'r')
        
        ng = max(HF1st(:));
        U_dis_in = Ucur(1:ng,7);

        
        %---- Generating Slices
        IntpMatrices2 = []; % set IntpMatrices empty only in the beginning.
        % Drawslices.m returns the IntpMatrices;
        % after this drawing the slices happens quickly
        
        xyslices2 = max(gtemp(:,3))*[0:.4:1];
        xzslices2 = max(gtemp(:,2))*[0];
        yzslices2 = max(gtemp(:,1))*[.5];
        
        figure(201)
        viewangle = [0 0];
        IntpMatrices2 = DrawSlices(gtemp,HF1st,NodeF1st,U_dis_in,xyslices2,xzslices2,yzslices2,'cylinder',max(gtemp(:,1))/50,201,IntpMatrices2,viewangle);
        colorbar
        title('with inclusion')

        
  % --- Difference
         
            
        Udiff = abs(U_dis_in - U_dis);          
        
        %---- Generating Slices
        IntpMatrices2 = []; % set IntpMatrices empty only in the beginning.
        % Drawslices.m returns the IntpMatrices;
        % after this drawing the slices happens quickly
        
        xyslices2 = max(gtemp(:,3))*[0:.4:1];
        xzslices2 = max(gtemp(:,2))*[0];
        yzslices2 = max(gtemp(:,1))*[.5];
        
        figure(101)
        viewangle = [0 0];
        IntpMatrices2 = DrawSlices(gtemp,HF1st,NodeF1st,Udiff,xyslices2,xzslices2,yzslices2,'cylinder',max(gtemp(:,1))/50,201,IntpMatrices2,viewangle);
        colorbar
        title('Difference')

        

