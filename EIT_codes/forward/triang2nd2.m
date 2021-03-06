function int=triang2nd2(g);

w=[1/40;1/15;1/40;1/15;1/40;1/15;9/40];
ip=[0 0; 1/2 0;1 0;1/2 1/2;0 1;0 1/2;1/3 1/3];

x=g(:,1);
y=g(:,2);
z=g(:,3);

int=0;
for ii=1:length(ip)
S=[2*ip(ii,1)^2+2*ip(ii,2)^2+4*ip(ii,1)*ip(ii,2)-3*ip(ii,1)-3*ip(ii,2)+1; ...
   2*ip(ii,1)^2-ip(ii,1); ...
   2*ip(ii,2)^2-ip(ii,2); ...
   -4*ip(ii,1)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,1); ...
   4*ip(ii,1)*ip(ii,2); ...
   -4*ip(ii,2)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,2)];

k=ip(ii,1);
n=ip(ii,2);
c=0;

aa=2*(2*x(1)+2*x(2)-4*x(4))*k+(-3*x(1)-x(2)+4*x(4))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*n;
bb=2*(2*y(1)+2*y(2)-4*y(4))*k+(-3*y(1)-y(2)+4*y(4))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*n;
cc=2*(2*z(1)+2*z(2)-4*z(4))*k+(-3*z(1)-z(2)+4*z(4))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*n;
AA=2*(2*x(1)+2*x(3)-4*x(6))*n+(-3*x(1)-x(3)+4*x(6))+(4*x(1)-4*x(4)+4*x(5)-4*x(6))*k;
BB=2*(2*y(1)+2*y(3)-4*y(6))*n+(-3*y(1)-y(3)+4*y(6))+(4*y(1)-4*y(4)+4*y(5)-4*y(6))*k;
CC=2*(2*z(1)+2*z(3)-4*z(6))*n+(-3*z(1)-z(3)+4*z(6))+(4*z(1)-4*z(4)+4*z(5)-4*z(6))*k;

ll=sqrt((bb*CC-BB*cc)^2+(cc*AA-CC*aa)^2+(aa*BB-AA*bb)^2);

int=int+w(ii)*ll*S*S';
 
end
