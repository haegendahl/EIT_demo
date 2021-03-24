function [Bubb] = GenerateAnisotropicRodInGrid(g,RANGE,WIDTH,POSITION,CONDUCTIVE_BLOB)

%
% 
% 

min_c = RANGE(1); 
max_c = RANGE(2);

width_x = WIDTH(1);
width_y = WIDTH(2);

xi = POSITION(1);
yi = POSITION(2);

Bubb = ones(length(g),1);
x=g(:,1);
y=g(:,2);

%keyboard

beta_x = width_x/6;
beta_y = width_y/6;

t = (1./sqrt(2*pi*beta_x^2)*exp(-((x-xi).^2)./(2*beta_x^2))).*...
    (1./sqrt(2*pi*beta_y^2)*exp(-((y-yi).^2)./(2*beta_y^2)));

BuB = 1/max(max(t))*t*.5;
Bubb = Bubb - BuB;


Bubb = (Bubb-min(min(Bubb)))/(max(max(Bubb))-min(min(Bubb)));
Bubb = min_c + (max_c-min_c)*Bubb;

if (strcmp(CONDUCTIVE_BLOB,'YES'))
  Bubb = min(min(Bubb))+max(max(Bubb)) - Bubb;
end


