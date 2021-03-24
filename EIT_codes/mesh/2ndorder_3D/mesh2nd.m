function [gn,Hn,elind,Eind] = mesh2nd(gn,Hn,elind,Eind,L)

%keyboard

gn1=[];
for i=1:length(Hn)
gn1=[gn1;gn(Hn(i,:),:)];
end

%%% lasketaan lisapisteet kulmapisteiden valiin %%%%%
gn2=[];
gn3=[];
for i=1:4:length(gn1)
 gn21(i,:)=(gn1(i,:)+gn1(i+1,:))/2;
 gn22(i,:)=(gn1(i+1,:)+gn1(i+2,:))/2;
 gn23(i,:)=(gn1(i,:)+gn1(i+2,:))/2;
 gn24(i,:)=(gn1(i,:)+gn1(i+3,:))/2;
 gn25(i,:)=(gn1(i+1,:)+gn1(i+3,:))/2;
 gn26(i,:)=(gn1(i+2,:)+gn1(i+3,:))/2;
 gn2=[gn2;gn21(i,:);gn22(i,:);gn23(i,:);gn24(i,:);gn25(i,:);gn26(i,:)];
 gn3=[gn3;gn1(i,:);gn21(i,:);gn1(i+1,:);gn22(i,:);gn1(i+2,:);gn23(i,:); ...
     gn1(i+3,:);gn24(i,:);gn25(i,:);gn26(i,:)];
end

%%%%%%% poistetaan samat pisteet %%%%%%%%%%%
gnlis=intersect(gn2,gn2,'rows');

gn=[gn;gnlis];

%%%%%%% topologia %%%%%%%%%%%%%%%%%%%%%%
%tic
%HHn=[];
%for jj=1:length(gn2)
%  [a,b,c]=intersect(gn2(jj,:),gn,'rows');
%   HHn=[HHn,c];
%end
%toc

tic
HHn=[];
gn31=gn2(:,1);
gn32=gn2(:,2);
gn33=gn2(:,3);
 for i=1:length(gn2)
   h=find(gn31(i)==gn(:,1));
   q=find(gn32(i)==gn(:,2));
   f=find(gn33(i)==gn(:,3));
    e=[];
    d=[];
    for j=1:length(h)
      e=[e;q(find(h(j)==q))];
      d=[d;f(find(h(j)==f))];
    end
     s=[];
     for k=1:length(e)
       s=[s;d(find(e(k)==d))];
     end
   HHn=[HHn;s];
   %i
end
toc

Hn1=[];
for i=1:6:length(HHn)
Hn1=[Hn1,HHn(i:i+5)];
end
Hn1=Hn1';
Hn=[Hn,Hn1];

%%%%%%%% elektrodit %%%%%%%%%%%%
elind1=[];
for ii=1:size(Eind,1)
 el=Eind(ii,:);
   for jj=1:length(el)
    h=ismember(Hn(el(jj),:),elind(:));
    if h(1)==1 & h(2)==1 & h(3)==1
     elind1=[elind1;Hn(el(jj),5:7)];
    end
    if h(1)==1 & h(2)==1 & h(4)==1
     elind1=[elind1;Hn(el(jj),[5 9 8])];
    end
    if h(2)==1 & h(3)==1 & h(4)==1
     elind1=[elind1;Hn(el(jj),[6 10 9])];
    end
    if h(1)==1 & h(3)==1 & h(4)==1
     elind1=[elind1;Hn(el(jj),[7 10 8])];
    end
  end
end

ee=reshape(elind1',3*size(Eind,2),L);
elind11=ee';

ell=[];
for jj=1:3:size(elind,2)
 ell=[ell,elind(:,jj:jj+2) elind11(:,jj:jj+2) ];
end
elind=ell;



