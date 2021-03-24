function [U,p,r]=ForwardSolution3dEl(Node,Element,T,f,sigma,z,MeasPattern,style);
% Same as ForwardSolution3d but computes ONLY the voltages on electrodes

% M. Vauhkonen 13.8.1999, University of Kuopio, Finland

p = 1; r =1;

L=max(size(z));
sN=max(size(Node));
[II1,C]=Current(L,sN,'adj');
C=C(:,1:L-1);% For the voltage reference 

if style=='comp'                    

%% CHANGED 12.01.2007
II1=sparse([zeros(L,sN),C,zeros(L,sN+L-1);zeros(L,2*sN+L-1),C]);
%II1=sparse([zeros(L,sN),C]);

if ~isempty(MeasPattern) 
 MeasPattern = [MeasPattern zeros(size(MeasPattern));zeros(size(MeasPattern)) MeasPattern];
 II1=MeasPattern'*II1;
 II1=II1'; 
else
 MeasPattern=eye(2*max(size(C)));
 II1=II1';
end

%% CHANGED 12.01.2007
II=sparse([[zeros(sN,size(T,2));C'*T];zeros(sN+L-1,size(T,2))]);
%II=[zeros(sN,size(T,2));C'*T];

% include angular frequency & permittivity of free space to sigma
w = 2*pi*f;
e0 = 8.85418781762e-12;
sigma = complex(real(sigma),w*e0*imag(sigma));
[A]=FemMatrix3d(Node,Element,sigma,z,C);

%% CHANGED 12.01.2007
A=[real(A),-imag(A);imag(A),real(A)];

UU=[A\II];%Voltages for the "measurement field" and for the current patterns.

U.Electrode=MeasPattern'*[C,zeros(size(C));zeros(size(C)),C]*[UU(sN+1:sN+L-1,:);UU(2*sN+L:size(UU,1),:)];%Voltages on the electrodes

elseif style=='real'
II1=sparse([zeros(L,sN),C]);
if ~isempty(MeasPattern) 
 II1=MeasPattern'*II1;
 II1=II1'; 
else
 MeasPattern=eye(max(size(C)));
 II1=II1';
end
II=[zeros(sN,size(T,2));C'*T];
[A]=FemMatrix3d(Node,Element,sigma,z,C);
%if nargin<8
 p=symamd(A);
 r(p)=1:max(size(p));
%end
%UU=A\[II1,II];%Voltages for the "measurement field" and for the current patterns.
R=chol(A(p,p));
UU=R\(R'\[II1(p,:),II(p,:)]);
UU=UU(r,:);

U.Electrode=MeasPattern'*C*UU(sN+1:size(A,1),:);%Voltages on the electrodes

end

