function [IntpMatrices] = DrawSlices(ginv,Hinv,Node,sigma,xyslices,xzslices,yzslices,gridtype,pixsize,figno,IntpMatrices,viewangle)

%keyboard

% Aku Seppänen 7.6.2007

nxyslices = length(xyslices);
nxzslices = length(xzslices);
nyzslices = length(yzslices);

totnrslices = nxyslices + nxzslices + nyzslices;

if(isempty(IntpMatrices))
  clear IntpMatrices
  for kk = 1:totnrslices
    IntpMatrices{kk} = sparse([]);
  end
end


minx = min(ginv(:,1));
maxx = max(ginv(:,1));
meanx = mean([minx maxx]);
miny = min(ginv(:,2));
maxy = max(ginv(:,2));
meany = mean([miny maxy]);
minz = min(ginv(:,3));
maxz = max(ginv(:,3));
meanz = mean([minz maxz]);
 
Nz = round((maxz-minz)/pixsize);
zstep = (maxz-minz)/Nz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Cylindrical                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(gridtype,'cylinder'))

    if(nxyslices > 0)

       maxr = maxx;
       rstep = round((maxr)/pixsize);
       rsteplngth = maxr/rstep;
       rr = (0:rsteplngth:maxr+10^-4)';
       thsteps = round(2*pi*rr/pixsize);
       thsteps(1) = 1;
       thsteplngths = 2*pi*rr(2:end)./thsteps(2:end);
       thsteplngths = [0; thsteplngths];
       %gtmpx = zeros(sum(thsteps),1)
       %gtmpy = zeros(sum(thsteps),1)
       gtmpx = [0];
       gtmpy = [0];
       for kk = 2:length(rr)
          th = [0:thsteplngths(kk):2*pi*rr(kk)-10^(-4)]'/rr(kk);
          rrr = rr(kk)*ones(size(th));
          [xtmp,ytmp] = pol2cart(th,rrr);
          gtmpx = [gtmpx; xtmp];
          gtmpy = [gtmpy; ytmp];
       end

       Htmp = delaunay(gtmpx,gtmpy);

       %figure,
       % trimesh(Htmp,gtmpx,gtmpy,ones(size(gtmpx,1),1),ones(size(gtmpx,1),1)),axis image
       %view(2)

       %gtmpzz = minz + (maxz-minz)*xydiv;

       for jj = 1:nxyslices

          gtmpz = xyslices(jj)*ones(size(gtmpx));
          VisNodes{jj} = [gtmpx gtmpy gtmpz];


          VisMatrixXY1 = IntpMatrices{jj};
          if(isempty(VisMatrixXY1))
	    [tmp VisMatrixXY1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma(:,end),VisNodes{jj},[]); 
          end

          VisElements{jj} = Htmp;
          IntpMatrices{jj} = VisMatrixXY1;

       end

    end


    if(nxzslices > 0)

       for jj = 1:nxzslices

	   minxjj = -sqrt(maxx^2 - (xzslices(jj))^2);
  	   maxxjj = -minxjj;

           Nxjj = round((maxxjj-minxjj)/pixsize);
           xstepjj = (maxxjj-minxjj)/Nxjj;
  	   [gtmpx gtmpz] = meshgrid(minxjj:xstepjj:maxxjj,minz:zstep:maxz);
           gtmpx = gtmpx(:);
           gtmpz = gtmpz(:);
  	   gtmpy = xzslices(jj)*ones(size(gtmpx));
           Htmp = delaunay(gtmpx,gtmpz);
           figure,
           trimesh(Htmp,gtmpx,gtmpz,ones(size(gtmpx,1),1),ones(size(gtmpx,1),1)),axis image
           view(2)
 
           VisNodes{nxyslices+jj} = [gtmpx gtmpy gtmpz];

           VisMatrixXZ1 = IntpMatrices{nxyslices+jj};
	   if(isempty(VisMatrixXZ1))
	     [tmp VisMatrixXZ1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma,VisNodes{nxyslices+jj},[]); 
           end

           VisElements{nxyslices+jj} = Htmp;
           IntpMatrices{nxyslices+jj} = VisMatrixXZ1;

       end

    end

    if(nyzslices > 0)

       for jj = 1:nyzslices

	   minyjj = -sqrt(maxy^2 - (yzslices(jj))^2);
  	   maxyjj = -minyjj;

           Nyjj = round((maxyjj-minyjj)/pixsize);
           ystepjj = (maxyjj-minyjj)/Nyjj;
  	   [gtmpy gtmpz] = meshgrid(minyjj:ystepjj:maxyjj,minz:zstep:maxz);
           gtmpy = gtmpy(:);
           gtmpz = gtmpz(:);
  	   gtmpx = yzslices(jj)*ones(size(gtmpy));
           Htmp = delaunay(gtmpy,gtmpz);
           %figure,
           %trimesh(Htmp,gtmpx,gtmpz,ones(size(gtmpx,1),1),ones(size(gtmpx,1),1)),axis image
           %view(2)
 
           VisNodes{nxyslices+nxzslices+jj} = [gtmpx gtmpy gtmpz];

           VisMatrixXZ1 = IntpMatrices{nxyslices+nxzslices+jj};
	   if(isempty(VisMatrixXZ1))
	     [tmp VisMatrixXZ1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma,VisNodes{nxyslices+nxzslices+jj},[]); 
           end

           VisElements{nxyslices+nxzslices+jj} = Htmp;
           IntpMatrices{nxyslices+nxzslices+jj} = VisMatrixXZ1;

       end

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Rectangular                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif(strcmp(gridtype,'rect'))

    if(nxyslices > 0)

	 minxjj = minx;
         maxxjj = maxx;
         Nxjj = round((maxxjj-minxjj)/pixsize);
         xstepjj = (maxxjj-minxjj)/Nxjj;
	 minyjj = miny;
  	 maxyjj = maxy;
         Nyjj = round((maxyjj-minyjj)/pixsize);
         ystepjj = (maxyjj-minyjj)/Nyjj;

  	 [gtmpx gtmpy] = meshgrid(minxjj:xstepjj:maxxjj,minyjj:ystepjj:maxyjj);
         gtmpx = gtmpx(:);
         gtmpy = gtmpy(:);

         Htmp = delaunay(gtmpx,gtmpy);

         for jj = 1:nxyslices

  	   gtmpz = xyslices(jj)*ones(size(gtmpx));
           VisNodes{jj} = [gtmpx gtmpy gtmpz];
           VisMatrixXZ1 = IntpMatrices{jj};
	   if(isempty(VisMatrixXZ1))
	     [tmp VisMatrixXZ1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma,VisNodes{jj},[]); 
           end

           VisElements{jj} = Htmp;
           IntpMatrices{jj} = VisMatrixXZ1;

         end

    end


    if(nxzslices > 0)

	 minxjj = minx;
         maxxjj = maxx;
         Nxjj = round((maxxjj-minxjj)/pixsize);
         xstepjj = (maxxjj-minxjj)/Nxjj;

  	 [gtmpx gtmpz] = meshgrid(minxjj:xstepjj:maxxjj,minz:zstep:maxz);
         gtmpx = gtmpx(:);
         gtmpz = gtmpz(:);
         Htmp = delaunay(gtmpx,gtmpz);

         for jj = 1:nxzslices

  	   gtmpy = xzslices(jj)*ones(size(gtmpx));
           VisNodes{nxyslices+jj} = [gtmpx gtmpy gtmpz];
           VisMatrixXZ1 = IntpMatrices{nxyslices+jj};
	   if(isempty(VisMatrixXZ1))
	     [tmp VisMatrixXZ1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma,VisNodes{nxyslices+jj},[]); 
           end

           VisElements{nxyslices+jj} = Htmp;
           IntpMatrices{nxyslices+jj} = VisMatrixXZ1;

       end

    end

    if(nyzslices > 0)

	   minyjj = miny;
  	   maxyjj = maxy;
           Nyjj = round((maxyjj-minyjj)/pixsize);
           ystepjj = (maxyjj-minyjj)/Nyjj;

  	   [gtmpy gtmpz] = meshgrid(minyjj:ystepjj:maxyjj,minz:zstep:maxz);
           gtmpy = gtmpy(:);
           gtmpz = gtmpz(:);
           Htmp = delaunay(gtmpy,gtmpz);
 
           for jj = 1:nyzslices

   	     gtmpx = yzslices(jj)*ones(size(gtmpy));
             VisNodes{nxyslices+nxzslices+jj} = [gtmpx gtmpy gtmpz];

             VisMatrixXZ1 = IntpMatrices{nxyslices+nxzslices+jj};
	     if(isempty(VisMatrixXZ1))
	       [tmp VisMatrixXZ1] = Interpolate2Newmesh3DNode(ginv,Hinv,Node,sigma,VisNodes{nxyslices+nxzslices+jj},[]); 
             end

             VisElements{nxyslices+nxzslices+jj} = Htmp;
             IntpMatrices{nxyslices+nxzslices+jj} = VisMatrixXZ1;

           end

      end

end

if(strcmp(gridtype,'cylinder') | strcmp(gridtype,'rect'))

  minval = min(min(sigma));
  maxval = max(max(sigma));

  figure(figno), clf, axis equal, hold on

  for ii = 1:length(VisNodes)
 
     H = VisElements{ii};
     g = VisNodes{ii};
     sig = IntpMatrices{ii}*sigma;
     fp = patch('faces',H,'vertices',g,'facevertexcdata',sig,'facecolor','interp','edgecolor','none');
     set(get(fp,'parent'),'CLim',[minval maxval])


  end
 %alpha(0.5)
  view(viewangle)


end
