% Draws boundary triangulation and electrodes
% in a cylindrical grid
% 
% Aku Sepp?nen 27.8.2002
%

figure(1), clf,hold on 
indBoundH=[];
for iH = 1:size(H,1)
     test=sqrt(g(H(iH,:),1).^2+g(H(iH,:),2).^2)>.99999999*R;
     if (sum(test)==3)
       indBoundH = [indBoundH;iH];
       indtest = find(test==1);
       indtest(4) = indtest(1);
       plot3(g(H(iH,indtest),1),g(H(iH,indtest),2),...
	     g(H(iH,indtest),3),'k')
     end
end
axis image
set(1,'Position',[349 204 633 645])
set(get(1,'child'),'CameraPosition',[178.211 -18.7307 96.2816])

% electrodes
for iel = 1:size(elind,1)
  for iH = 1:size(elind,2)/3
     gdraw = g(elind(iel,(iH-1)*3+1:iH*3),:);
     patch(gdraw(:,1),gdraw(:,2),gdraw(:,3),'k')
  end
end





