function [Fnew,Urefel,Uref] = costfunTVnoneg(theta,cost_params)
   
Node     = cost_params.Node;
Element  = cost_params.Element;
I        = cost_params.I;
P1st     = cost_params.P1st;
z        = cost_params.z;
MeasPatt = cost_params.MeasPatt;
Uel      = cost_params.Uel;
Ln       = cost_params.Ln;
R        = cost_params.R;
beta     = cost_params.beta;
alpha    = cost_params.alpha;
a        = cost_params.a;


Uref = ForwardSolution3d2ndElectrode_fix(Node,Element,I,P1st*theta,z,MeasPatt,'real');
Urefel = Uref(:);


Rt = R*theta;
Rx  = Rt(1:end/3,:);
Ry = Rt(end/3+1:2*end/3,:);
Rz  = Rt(1+2*(end/3):3*end/3,:);
t = sum(sqrt( Rx.^2 + Ry.^2 + Rz.^2 + beta));

q = -sum(a.*log(theta));

Fnew = .5*norm(Ln*(Uel-Urefel))^2 + alpha*t + q;



