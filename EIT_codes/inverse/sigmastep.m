function [beta,I,J] = sigmastep(sigma,d,epsilon_s);
% function assumes sigma as a positive double

% epsilon_s = greatest lower bound for conductivies (if direction of the conductitivity is
% towards the zero, the conductivity won't be updated after this bound is reached)


ns = length(sigma(:));
if nargin<3, epsilon_s = 1e-3; end

epsilon = epsilon_s*.0025;     % lowest possible conductivity to be achieved with a single step
if min(sigma)<epsilon
  epsilon = min(sigma);
end

beta = ((epsilon + 0.01*epsilon) - sigma)./d;  % step lengths
I = find((sigma>epsilon_s) & (beta>0));

if isempty(I)  % all parameters move towards the +\infty
  bmax = 1;
else           % some move towards the zero (-\infty)
  bmax = min(beta(I));
end

beta = (epsilon_s - sigma)./d;
J = find((beta<0) & (sigma<=epsilon_s));
I = setdiff([1:ns]',J);


beta(I) = bmax;
beta(J) = 0;
