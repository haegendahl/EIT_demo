function [in]=isinside(p1,p2,p3,p4,p);
%function [in]=isinside(p1,p2,p3,p4,p);
%Chechks if point is inside tetrahedra
% if inside in=1 else in=0
%ponts p1 to p4 are vertice points of tetrahedra
% and point p is the point under inspection
% tason yht. Ax+By+Cz+D==0
% d is distance of point p from taso
in=0;
in1=0;in2=0;in3=0;in4=0;
x=[p1(1) p2(1) p3(1)]';
y=[p1(2) p2(2) p3(2)]';
z=[p1(3) p2(3) p3(3)]';

A=((y(2)-y(1))*(z(3)-z(1))-(y(3)-y(1))*(z(2)-z(1)));
B=-((x(2)-x(1))*(z(3)-z(1))-(x(3)-x(1))*(z(2)-z(1)));
C=((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)));
D=(-x(1)*A-y(1)*B-z(1)*C);

dp4=(A*p4(1)+B*p4(2)+C*p4(3)+D)/sqrt(A*A+B*B+C*C);
%from the sign of dp4 we get the sign for dp
%it must have same sign if it's inside tetrahedra
dp=(A*p(1)+B*p(2)+C*p(3)+D)/sqrt(A*A+B*B+C*C);
if (dp*dp4 >= 0) & (abs(dp)<=abs(dp4)) 
 in1=1; else in1=0;end; 

if in1==1
in2=0;
x=[ p2(1) p3(1) p4(1) ]';
y=[ p2(2) p3(2) p4(2)]';
z=[ p2(3) p3(3) p4(3)]';

A=((y(2)-y(1))*(z(3)-z(1))-(y(3)-y(1))*(z(2)-z(1)));
B=-((x(2)-x(1))*(z(3)-z(1))-(x(3)-x(1))*(z(2)-z(1)));
C=((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)));
D=(-x(1)*A-y(1)*B-z(1)*C);

dp1=(A*p1(1)+B*p1(2)+C*p1(3)+D)/sqrt(A*A+B*B+C*C);
%from the sign of dp4 we get the sign for dp
%it must have same sign if it's inside tetrahedra
dp=(A*p(1)+B*p(2)+C*p(3)+D)/sqrt(A*A+B*B+C*C);
if (dp*dp1 >= 0) & (abs(dp)<=abs(dp1)) 
 in2=1; else in2=0;end;
end
 
if (in1==1) & (in2==1)
in3=0;
x=[ p1(1) p3(1) p4(1) ]';
y=[ p1(2) p3(2) p4(2)]';
z=[ p1(3) p3(3) p4(3)]';

A=((y(2)-y(1))*(z(3)-z(1))-(y(3)-y(1))*(z(2)-z(1)));
B=-((x(2)-x(1))*(z(3)-z(1))-(x(3)-x(1))*(z(2)-z(1)));
C=((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)));
D=(-x(1)*A-y(1)*B-z(1)*C);
dp2=(A*p2(1)+B*p2(2)+C*p2(3)+D)/sqrt(A*A+B*B+C*C);
dp=(A*p(1)+B*p(2)+C*p(3)+D)/sqrt(A*A+B*B+C*C);
if (dp*dp2 >= 0) & (abs(dp)<=abs(dp2)) 
 in3=1; else in3=0;end;
end

if (in1==1) & (in2==1) & (in3==1)
in4=0;
x=[ p1(1) p2(1) p4(1) ]';
y=[ p1(2) p2(2) p4(2)]';
z=[ p1(3) p2(3) p4(3)]';

A=((y(2)-y(1))*(z(3)-z(1))-(y(3)-y(1))*(z(2)-z(1)));
B=-((x(2)-x(1))*(z(3)-z(1))-(x(3)-x(1))*(z(2)-z(1)));
C=((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)));
D=(-x(1)*A-y(1)*B-z(1)*C);
dp3=(A*p3(1)+B*p3(2)+C*p3(3)+D)/sqrt(A*A+B*B+C*C);
dp=(A*p(1)+B*p(2)+C*p(3)+D)/sqrt(A*A+B*B+C*C);
if (dp*dp3 >= 0) & (abs(dp)<=abs(dp3)) 
 in4=1; else in4=0;end;
end

[in1 in2 in3 in4];

if (in1==1) & (in2==1) & (in3==1) & (in4==1)
	in=1;
else
	in=0;
end


