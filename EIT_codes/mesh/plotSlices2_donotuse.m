function plotSlices2(sigma,h,g,nz,z)
% simple slice plot
%disp('lol')
ng = size(g,1);
ns = size(sigma,1);

if mod(ns,ng)
  error('lol!')
  return
end
 
nn = ns/ng;

hold on
alphaslice = ones(nn,1);
%keyboard
for ii=nz
  ss = sigma((ii-1)*ng+1:ii*ng);
  patch('faces',h,'vertices',[g z(ii)*ones(size(g,1),1)],'facevertexcdata',ss,'facecolor','interp','edgecolor','none','facealpha',alphaslice(ii));
end
axis equal, view(3)
