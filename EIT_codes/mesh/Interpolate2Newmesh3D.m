function [f_newgrid,INTPMAT,Emap] = Interpolate2Newmesh3D(g,H,f,pts,INTPMAT,Emap)

%keyboard

disp('Fast sparse!')
if(isempty(INTPMAT))

 
   

   centers = zeros(size(H,1),3);
   for jj = 1:size(H,1)
        centers(jj,:) = mean(g(H(jj,:),:),1);
   end;

   invX = cell(size(H,1),1);
   for tin = 1:size(H,1);

        pp = H(tin,:);
        X = [g(pp,:),ones(4,1)];
        invX{tin} = inv(X);
 
   end;


   np = size(pts,1);
   Ic = zeros(np,4);
   Iv = zeros(np,4);

   createEmap = 0;
   if isempty(Emap), Emap = cell(np,1); createEmap=1; end
   for k = 1:np

      %%%% determine the tedrahedron %%%% 

      x = pts(k,1); y = pts(k,2); z = pts(k,3);
      if createEmap
	%tin = FindElement(g,H,[0 0 0; pts(k,:)]);
        tin = tsearchn(g,H,pts(k,:));
	Emap{k} = tin;
      else
	tin = Emap{k};
      end
  
      if(isempty(tin))
        disp('fuck')
        dst = (centers(:,1)-x).^2+(centers(:,2)-y).^2+(centers(:,3)-z).^2;
        tin = find(dst==min(dst));
        %tin = tin(1); %jos yhtä lähellä kahta keskipistettä
      end;
      tin = tin(1); %jos yhtä lähellä kahta keskipistettä


      for gin = 1:4
           I = zeros(4,1);
           I(gin)=1; 
           Phi(gin) = [x y z 1]*invX{tin}*I;
      end;

      Ic(k,:) = H(tin,:);
      Iv(k,:) = Phi;      

   end;
   INTPMAT = sparse(repmat([1:np]',1,4),Ic,Iv,np,size(f,1));
   % INTPMAT = sparse(size(pts,1),size(f,1));
   
end;

f_newgrid = INTPMAT*f;










