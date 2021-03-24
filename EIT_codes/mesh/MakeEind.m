function  [Eind,elind]=MakeEind(gn,Hn,E,ntas,N,lg,lindb);

% M. Vauhkonen 2.12.1999, University of Kuopio, Dept. of Applied Physics
% Modified by Aku Seppänen 27.8.2002

%keyboard

gind=[];
elind=[];
kk=1;
Eind=[];
ek1=[];
ek2=[];


for ii=1:max(size(N))
  if N(ii)~=0
   n=ii-1;
   ek1=[3*E-2+n*length(Hn)/(ntas-1)];
   ek2=[3*E-1+n*length(Hn)/(ntas-1)];
   gind=[gind;[n*lg+1:n*lg+lindb]];
   Eind=[Eind;[ek1,ek2]];
   gind=[gind;(n+1)*lg+1:(n+1)*lg+lindb];
  end 
end

 Eind1=Eind';
 Eind1=Eind1(:);
 

A = Hn(Eind1,:);
a=ismember(A,gind);

elind=[];
for ii=1:length(a)
  [b,c] = find(a(ii,:)==1);
  elind = [elind;A(ii,c)];
end  
%elind = (reshape(elind',3*size(Eind,2),size(Eind,1)))';
%size(elind)

%%%%% Do this ONLY for the FEM calculations...
%L=size(find(N),1)*max(N); %Number of electrodes.
L=size(E,1)*sum(N); %Number of electrodes.
elind=elind';
no=3; %number of nodes in triangle
elind=reshape(elind,size(elind,2)/L*no,L);
%keyboard
elind=elind';

% If there are electrodes wider than one layer (N=[... 0 1 1 0 ..])
N_ones = find(N);
N_zeros = find(N==0);
tmp = find(N_zeros>N_ones(1));
w_el_layer = N_zeros(tmp(1)) - N_ones(1);  % width of the electrode layers
N_el_layer = length(N_ones)/w_el_layer;    % number of the electrode layers
Nel = L/w_el_layer;                        % actual number of electrodes
elind_new = zeros(N_el_layer,size(elind,2)*w_el_layer);
Eind_new = zeros(N_el_layer,size(Eind,2)*w_el_layer);
for iel = 1:Nel
   ilayer = ceil(iel/size(E,1));
   iel_l = iel-(ilayer-1)*size(E,1);
   iiel = (ilayer-1)*size(E,1)*w_el_layer + iel_l:size(E,1):...
            (ilayer-1)*size(E,1)*w_el_layer + iel_l + size(E,1)*(w_el_layer-1);
   elind_new_iel = (elind(iiel,:))';
   elind_new(iel,:) = (elind_new_iel(:))';
   Eind_new_iel = (Eind(iiel,:))';
   Eind_new(iel,:) = (Eind_new_iel(:))';
end

elind = elind_new;
Eind = Eind_new;














































































































