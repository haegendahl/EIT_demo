function [fh,fp]=PlotSolScaledImage(g,H,V,minval,maxval,N);

% Aku Seppänen 24.11.1999.

%keyboard

[VN,VM]=size(V);
[gN,gM]=size(g);
[HN,HM]=size(H);

% New rectangular grid
mn_x = min(g(:,1)); mx_x = max(g(:,1));
mn_y = min(g(:,2)); mx_y = max(g(:,2));
dx = (mx_x-mn_x)/N;
x = (mn_x:dx:mx_x)';y=(mn_y:dx:mx_y)';
[xx,yy] = meshgrid(x,y);
xx=xx(:); yy=yy(:);

% coefficients of basis functions
V_pol_c = zeros(length(H),9);

for tin = 1:length(H)
     pp=H(tin,:);
     x_j = g(pp(1),1); y_j = g(pp(1),2);
     x_l = g(pp(2),1); y_l = g(pp(2),2);
     x_k = g(pp(3),1); y_k = g(pp(3),2);
     alpha_j = V(pp(1));
     alpha_l = V(pp(2));
     alpha_k = V(pp(3));

     if x_k==x_l
       b_j=0;
       a_j=alpha_j/(x_j-x_k);
       c_j=-a_j*x_k;
     else
       b_j=alpha_j/((y_l-y_k)/(x_k-x_l)*(x_j-x_k)+y_j-y_k);
       a_j=b_j*(y_l-y_k)/(x_k-x_l);
       c_j=-b_j*((y_l-y_k)/(x_k-x_l)*x_k+y_k);
     end;

     if x_l==x_j
       b_k=0;
       a_k=alpha_k/(x_k-x_l);
       c_k=-a_k*x_l;
     else
       b_k=alpha_k/((y_j-y_l)/(x_l-x_j)*(x_k-x_l)+y_k-y_l);
       a_k=b_k*(y_j-y_l)/(x_l-x_j);
       c_k=-b_k*((y_j-y_l)/(x_l-x_j)*x_l+y_l);
     end;
    
     if x_j==x_k
       b_l=0;
       a_l=alpha_l/(x_l-x_j);
       c_l=-a_l*x_j;
     else
       b_l=alpha_l/((y_k-y_j)/(x_j-x_k)*(x_l-x_j)+y_l-y_j);
       a_l=b_l*(y_k-y_j)/(x_j-x_k);
       c_l=-b_l*((y_k-y_j)/(x_j-x_k)*x_j+y_j);
     end;

     V_pol_c(tin,:) = [a_j b_j c_j a_l b_l c_l a_k b_k c_k]; 

  end;


%%% Interpolate the values of function V in points x,y %%%

V_p = zeros(length(xx),1);
for k = 1:length(xx)
   p = [xx(k),yy(k)];
   h = [p(1),p(2),1,p(1),p(2),1,p(1),p(2),1];
   t_in = tsearch(g(:,1),g(:,2),H,p(1),p(2));
   if (find(t_in>0)>0)
     Pol=(V_pol_c(t_in,:))';
     V_p(k) = h*Pol;
   else
     V_p(k) = NaN;
   end;
end;

%keyboard

colormap(gray);
%map = colormap;
%map(1,:) = [1 1 1]; %white
%colormap(map);
dv = (maxval - minval)/(length(map));
imagesc(flipud(reshape(V_p,length(y),length(x))),[minval-dv,maxval]),axis image,axis off

 

