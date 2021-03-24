function plotSlices(sigma,h,g,z)
% simple slice plot
%disp('lol')
ng = size(g,1);
ns = size(sigma,1);

if mod(ns,ng)
  error('lol!')
  return
end

nn = ns/ng

pos = 1;

hold on
alphaslice = [1 1 1 1 1 1 1 1];
for ii=1:7
  ss = sigma(pos:pos+ng-1);
  patch('faces',h,'vertices',[g z(ii)*ones(size(g,1),1)],'facevertexcdata',ss,'facecolor','interp','edgecolor','none','facealpha',alphaslice(ii));
  pos = pos + ng;
end
axis equal, view(3)
