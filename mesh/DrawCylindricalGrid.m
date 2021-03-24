function [] = DrawCylindricalGrid(g,H,elind)
% Draws boundary triangulation and electrodes
% in a cylindrical grid
% 
% Aku Seppänen 27.8.2002
% Aku Seppänen 21.2.2008

%keyboard

figure(101), clf,hold on 
indBoundH=[];
for iH = 1:size(H,1)
     test = sqrt(g(H(iH,:),1).^2+g(H(iH,:),2).^2)>.99999999*max(g(:,1));
     test2 = g(H(iH,:),3)-10^(-8) < min(g(:,3));
     test3 = g(H(iH,:),3)+10^(-8) > max(g(:,3));

     if (sum(test)==3)
       indBoundH = [indBoundH;iH];
       indtest = find(test==1);
       indtest(4) = indtest(1);
       plot3(g(H(iH,indtest),1),g(H(iH,indtest),2),...
	     g(H(iH,indtest),3),'k')
     end

     if (sum(test2)==3)
       indBoundH = [indBoundH;iH];
       indtest = find(test2==1);
       indtest(4) = indtest(1);
       plot3(g(H(iH,indtest),1),g(H(iH,indtest),2),...
	     g(H(iH,indtest),3),'y')
     end

     if (sum(test3)==3)
       indBoundH = [indBoundH;iH];
       indtest = find(test3==1);
       indtest(4) = indtest(1);
       plot3(g(H(iH,indtest),1),g(H(iH,indtest),2),...
	     g(H(iH,indtest),3),'k')
     end

end



axis image
set(101,'Position',[349 204 633 645])
set(get(101,'child'),'CameraPosition',[178.211 -18.7307 96.2816])

% electrodes
for iel = 1:size(elind,1)
  for iH = 1:size(elind,2)/3
     gdraw = g(elind(iel,(iH-1)*3+1:iH*3),:);
     patch(gdraw(:,1),gdraw(:,2),gdraw(:,3),'k')
  end
end





