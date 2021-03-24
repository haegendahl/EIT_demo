function DrawCircle(x,y,z,r)

ang = [linspace(0,2*pi,100)]';

xp = r*cos(ang);

yp = r*sin(ang);

z = z*ones(size(ang,1),1);

plot3(x+xp,y+yp,z,'k','Color','black','LineWidth',1);