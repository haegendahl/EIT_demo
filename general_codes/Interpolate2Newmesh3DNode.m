function [f_newgrid,INTPMAT] = Interpolate2Newmesh3DNode(g,H,Node,f,pts,INTPMAT)

%keyboard

disp('Fast sparse!')
if(isempty(INTPMAT))

% $$$   centers = zeros(size(H,1),3);
% $$$   no_fcuk = 1; 
% $$$   
   invX = cell(size(H,1),1);
   for tin = 1:size(H,1);

        pp = H(tin,:);
        X = [g(pp,:),ones(4,1)];
        invX{tin} = inv(X);
 
   end;


   np = size(pts,1);
   Ic = zeros(np,4);
   Iv = zeros(np,4);

   [Element,minnode] = FindElementNode(g,H,Node,pts);
   
   for k = 1:np

      %%%% determine the tedrahedron %%%% 

      x = pts(k,1); y = pts(k,2); z = pts(k,3);
      tin = Element(k);
      
      Phi = zeros(1,4);
      if tin       
        iXt = invX{tin};
        for gin = 1:4
          Phi(gin) = [x y z 1]*iXt(:,gin);
        end;
        Ic(k,:) = H(tin,:);
      else  
        Ic(k,:) = [minnode(k) 1 1 1];
        Phi(1) = 1;
      end
            
      Iv(k,:) = Phi;      

   end;
   INTPMAT = sparse(repmat([1:np]',1,4),Ic,Iv,np,size(f,1));
   % INTPMAT = sparse(size(pts,1),size(f,1));
   
end;

f_newgrid = INTPMAT*f;
% $$$ 
% $$$ 
% $$$        if no_fcuk
% $$$           for jj = 1:size(H,1)
% $$$             centers(jj,:) = mean(g(H(jj,:),:),1);
% $$$           end;
% $$$           no_fcuk = 0;
% $$$         end        
% $$$         
% $$$         dst = (centers(:,1)-x).^2+(centers(:,2)-y).^2+(centers(:,3)-z).^2;
% $$$         [dst tin] = min(dst);
% $$$         
% $$$       end;
% $$$       tin = tin(1); %jos yhtä lähellä kahta keskipistettä
% $$$ 
% $$$ 






