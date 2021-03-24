function [I,T]=current_eit(L,lg,style,gn);

% Function [I,T]=current_eit(L,lg,style) calculates 
% different current patterns. style can be 'tri', 
% 'adj' or 'ref', L is number of the electrodes
% and lg is the number of the nodes. 

% M. Vauhkonen 6.9.1994. Modified 13.9.1994 for NOSER.

if style=='tri'
II=zeros(L,L-1);
l=(1:L)';
th=2*pi*l/L;
 for k=1:L/2,
  II(:,k)=cos(k*th);
 end
 for k=L/2+1:L-1
  II(:,k)=sin((k-L/2)*th);
 end
I=[sparse(zeros(lg,L-1));II];
T=II;
end

if style=='adj'
%II=zeros(L,L-1);
II=zeros(L,L);
l=(1:L)';
II1=diag(ones(L,1));
II2=diag(ones(L-1,1),-1);
II=II1(:,1:L-1)-II2(:,1:L-1);
II=[II,[-1;zeros(L-2,1);1]]; % lissee
%I=[zeros(lg,L-1);IIL];
I=[zeros(lg,L);II];
T=II;
end

if style=='ref'
%II=zeros(L,L-1);
II=zeros(L,L-1);
II1=ones(1,L-1);
II2=-1*diag(ones(L-1,1));
II=[II1;II2];
I=[zeros(lg,L-1);II];
T=II;
end

if style=='a3r'
%a=0.000264/sqrt(2);
%a=1;
a=0.001;
f=16;
I=zeros(f,f);
l=(1:f)';
II1=a*diag(ones(f,1));
II2=a*diag(ones(f-1,1),-1);
II=II1(:,1:f-1)-II2(:,1:f-1);
II=[II,[-a;zeros(f-2,1);a]]; % lissee
tt=zeros(f,f);
tt1=[II,tt,tt];
tt2=[tt,II,tt];
tt3=[tt,tt,II];
tot=[tt1;tt2;tt3];
I=[zeros(lg,L);tot]; %length(gn)=lg
T=tot;
end

if style=='a5r'
%a=0.000264/sqrt(2);
a=0.001;
%a=1;
f=16;
I=zeros(f,f);
l=(1:f)';
II1=a*diag(ones(f,1));
II2=a*diag(ones(f-1,1),-1);
II=II1(:,1:f-1)-II2(:,1:f-1);
II=[II,[-a;zeros(f-2,1);a]]; % lissee
tt=zeros(f,f);
tt1=[II,tt,tt,tt,tt];
tt2=[tt,II,tt,tt,tt];
tt3=[tt,tt,II,tt,tt];
tt4=[tt,tt,tt,II,tt];
tt5=[tt,tt,tt,tt,II];
tot=[tt1;tt2;tt3;tt4;tt5];
I=[zeros(lg,L);tot]; %length(gn)=lg
T=tot;
end

if style=='a6r'
%a=0.000236/sqrt(2);
a=0.000264/sqrt(2);
%a=1;
f=16;
II=zeros(f,f);
l=(1:f)';
II1=a*diag(ones(f,1));
II2=a*diag(ones(f-1,1),-1);
II=II1(:,1:f-1)-II2(:,1:f-1);
II=[II,[-a;zeros(f-2,1);a]]; % lissee
tt=zeros(f,f);
tt1=[II,tt,tt,tt,tt,tt];
tt2=[tt,II,tt,tt,tt,tt];
tt3=[tt,tt,II,tt,tt,tt];
tt4=[tt,tt,tt,II,tt,tt];
tt5=[tt,tt,tt,tt,II,tt];
tt6=[tt,tt,tt,tt,tt,II];
tot=[tt1;tt2;tt3;tt4;tt5;tt6];
I=[zeros(lg,L);tot]; %length(gn)=lg
 T=tot;
end


if style=='cee'
I=[];
C1=diag(ones(L,1));
C2=-diag(ones(L,1),-1);
T=C1+C2(1:L,1:L);
end

%if style=='cee2'
%T=zeros(48,96);
%C1=diag(ones(L/3,1));
%C2=-diag(ones(L/3,1),-2);
%TT=C1+C2(1:L/3,1:L/3);
%T(1:16,1:16)=TT;
%T(1:16,17:32)=TT;
%T(17:32,33:48)=TT;
%T(17:32,49:64)=TT;
%T(33:48,65:80)=TT;
%T(33:48,81:96)=TT;
%end

%if style=='6rjt'
%a=[1 0 -1]';
%C1=zeros(16,16);
%for i=1:8
%C1(i*2-1:i*2+1,i)=a;
%end
%C1(1,:)=C1(1,:)+C1(17,:);
%C1(17,:)=[];
%aa=zeros(3*16,6*16);
%for i=0:5
%aa(i*16+1:(i+1)*16,i*16+1:(i+1)*16)=C1;
%end
%I=[zeros(lg,L);aa]; %length(gn)=lg
%T=aa;
%end































