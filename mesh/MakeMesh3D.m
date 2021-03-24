function [Eind,elind,gn,Hn,g,H]=makegrid3new(N,r,NN,h,el,zs);

% Function [Eind,elind,gn,Hn]=makegrid3new(el,N,r,NN,h);
% calculates a 3D cylinderical grid, give parameters
% zs = separation of the layers in z-direction
% N = electrode layer/empty layer definitions
% e.g. N=[0,16,0,16,0,16,0]
% r = Radia of the circles, e.g. r=[1,0.7,0.35,0];
% NN = Number of nodes on each radius e.g. NN=[16,8,4,1];
% h = The height of the mesh
% el =electrode definition

% M. Vauhkonen 28.10. 1999, University of Kuopio, Dept. of Applied Physics, Finland


ntas=max(size(N))+1;  % ntas = The number of node layers in z-direction
Ne=max(N);

[g,gp,H,E]=cirgrid_eit(r,NN,el);
kk1=find((H(:,2)-H(:,1))<0);
kk2=find((H(:,3)-H(:,2))<0);
kk=[kk1;kk2];
lg=max(size(g));

Hn=[];lg=length(g);
for ii=1:length(H)
 ind=(H(ii,:));
 gpp=gp(ind,:);
 if length(find(gpp(:,1)==max(gpp(:,1))))==2 % Two nodes in the outer layer.
  ind1=find(gpp(:,1)==max(gpp(:,1)))';
  mind=ind(ind1(1));maind=ind(ind1(2));
  
 if any(ii==kk)
  ind=sort(H(ii,:));
  mind=ind(2);maind=ind(1); 
 end
  Hn1=[mind,mind+lg,maind+lg,max(ind); ...
      mind,maind,maind+lg,max(ind); ...
      mind+lg,maind+lg,max(ind)+lg,max(ind)]; 
  Hn=[Hn;Hn1];
 else
  ind1=find(gpp(:,1)==min(gpp(:,1)))';
  mind=ind(ind1(1));maind=ind(ind1(2)); 
 if any(ii==kk)
  ind=sort(H(ii,:));
  mind=ind(3);maind=ind(2);
 end
  Hn1=[mind,maind,maind+lg,min(ind)+lg; ...
      mind,maind,min(ind)+lg,min(ind); ...
      mind,mind+lg,maind+lg,min(ind)+lg];
  Hn=[Hn;Hn1];
 end
end 

%%%% Calculate all the nodes...
%z=[0:h/(ntas-1):h];
z=[0;cumsum(zs)];
gn=[];
 for ii=1:length(z)
  gn=[gn;[g,z(ii)*ones(length(g),1)]];
 end

%%%%% Calculate all the elements...

Hnl1=Hn;
for ii=1:ntas-2
 Hn=[Hn;Hnl1+ii*(lg)];
end

%indb=findboundary(g,H);
%lindb=size(indb);
[Eind,elind]=MakeEind(gn,Hn,E,ntas,N,lg,NN(1));




























































