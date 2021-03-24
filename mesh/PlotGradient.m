function PlotGradient(THETA, R, g,H,figno)


Rthta = R*THETA;

norm_thta = sqrt(Rthta(1:end/2,:).^2 +Rthta(1+end/2:end,:).^2 );


figure(figno), clf
axis off
axis equal
markgrid(norm_thta,g,H);

