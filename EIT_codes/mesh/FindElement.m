function [dipind] = FindElement(g,H,d);
%function [dipind]=dipelem3(g,H,d);
% TEMPORARY VERSION NOT QUITE CORRECT
% finds element on which the dipole is.
% d has form [Dx,Dy,Dz ; rox,roy,roz]
% dipole separation is infit. small so
% only location (rox,roy,roz) is checked 
% calls function isinside which checks if
% the point is inside tetrahedra

x=d(2,1); y=d(2,2); z=d(2,3);

Hind=[];
gmin=sqrt((g(:,1)-x).^2 + (g(:,2)-y).^2 + (g(:,3)-z).^2);
gmin=find(gmin==min(gmin));
gmin;
gind=gmin(1);
%haetaan pisteessa kiinni olevat tetraedrit

Hind=find( H(:,1)==gind | H(:,2)==gind | H(:,3)==gind | H(:,4)==gind);

dipind=[];
%tarkistetaan onko piste elementin sisalla
di=0;yes=0;
for ii=1:length(Hind)
yes=isinside(g(H(Hind(ii),1),:), g(H(Hind(ii),2),:) ,...
 g(H(Hind(ii),3),:),g(H(Hind(ii),4),:) ,[x,y,z]);
if (yes==1) di=di+1;
dipind=Hind(ii);end;
yes=0;
end;



