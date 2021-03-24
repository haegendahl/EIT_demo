function [I,T]=Current(L,lg,style,rms,numpat);

% Function [I,T]=Current(L,lg,style) calculates 
% different current patterns. style can be 'tri' (trigonometric), 
% 'adj' (adjacent), 'opp' (opposite) or 'ref' (one reference). 
%
% L = the number of electrodes
% lg = the number of nodes. 
% rms = the root mean square value of the injected current. Default=1.
% numpat= the number of patterns

% M. Vauhkonen 6.9.1994. Modified 13.9.1994 for NOSER.

I = 1; T = 1;


if nargin < 4 & style == 'tri'
 rms=1;
 numpat=L-1;
end

if nargin < 4 & style ~= 'tri'
 rms=1;
 numpat=L;
end

if nargin == 4 & style == 'tri'
 numpat=L-1;
end

if nargin == 4 & style ~= 'tri'
 numpat=L;
end


switch style
case 'tri'
if numpat > 15,
  str=['Too many current patterns. Automatically set to ', num2str(L-1),'.' ];
  disp(str)
  numpat =15;
 end
II=zeros(L,numpat);
l=(1:L)';
th=2*pi*l/L;
 for k=1:(numpat+1)/2,
  II(:,k)=rms*cos(k*th);
 end
 for k=(numpat+1)/2+1:numpat
  II(:,k)=rms*sin((k-L/2)*th);
 end
%nor=sqrt(diag(II'*II)*ones(1,L));
%II=II./nor';
I=[zeros(lg,numpat);II];
T=II;

case 'adj'
II=zeros(L,numpat);
l=(1:L)';
II1=diag(ones(L,1));
II2=diag(ones(L-1,1),-1);
if numpat < L
 II=II1(:,1:numpat)-II2(:,1:numpat);
else
 II=II1(:,1:numpat-1)-II2(:,1:numpat-1);
 II=[II,[-1;zeros(L-2,1);1]];% Optional!!
end
I=[zeros(lg,numpat);rms*II];
T=II;

case 'ref'
 if numpat > 15, 
  str=['Too many current patterns. Automatically set to ', num2str(L-1),'.' ];
  disp(str)
  numpat =15;
 end
II=zeros(L,numpat);
II1=diag(ones(L-1,1),-1);
II(1,:)=ones(1,numpat);
T=II-II1(:,1:numpat);
I=[zeros(lg,numpat);rms*T];


case 'opp'
if numpat > L/2, 
str=['Too many current patterns. Automatically set to ', num2str(L/2),'.' ];
  disp(str)
numpat=L/2;end
II=zeros(L,numpat);
l=(1:L)';
II1=diag(ones(L,1));
II2=diag(ones(L,1),-8);
II=II1(:,1:numpat)-II2(1:L,1:numpat);
%II=[II(:,1:ceil(numpat/2)),-II(:,1:floor(numpat/2))];% Optional!!
I=[zeros(lg,numpat);rms*II];
T=II;

case 'cee'
 %Note, L is here the total number of electrodes
 I=[];
 C1=diag(ones(L,1));
 C2=-diag(ones(L,1),-1);
 T=C1+C2(1:L,1:L);

case 'a2r'
 %Note, L is here the number of electrodes on each layer
 I=zeros(L,L);
 l=(1:L)';
 II1=rms*diag(ones(L,1));
 II2=rms*diag(ones(L-1,1),-1);
 II=II1(:,1:L-1)-II2(:,1:L-1);
 II=[II,[-rms;zeros(L-2,1);rms]]; % lissee
 tt=zeros(L,L);
 tt1=[II,tt];
 tt2=[tt,II];
 tot=[tt1;tt2];
 I=[zeros(lg,size(tot,2));tot]; %length(gn)=lg
 T=tot;
end










