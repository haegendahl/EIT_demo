function []=markgrid(rho,g,H);
% Need: N = size of the grid, g = grid point matrix, 
% H = index matrix of points of g and Igbord,Igint = grid point indexes 
% of the border and interior of the circle, respectively.
% C(size(H,1),1) is color for each triangle

% J. Kaipio, 11.4.1994. Modified by J. Virkkala 19.9.1994

% Callback strings

% Plot the existing basis
t=[0:.1:2.1*pi]';
[nH,mH]=size(H);

%plot(sin(t),cos(t),'g')
%axis([-1.2 1.2 -1.2 1.2])
%axis(axis)
%hold on
% plot(g(:,1),g(:,2),'w.');
% hbordpts=plot(g(Igbord,1),g(Igbord,2),'wo');

%colormap('jet')
%axis('off'),axis('square')
set(gcf,'defaultpatchedgecolor','none');	% CINVERT

%hdl=[];
for ii=1:nH % nH
  Hii=g(H(ii,:),:);
%  hdl=[hdl 
   patch(Hii(:,1),Hii(:,2),rho(ii));		% NEW, PATCH HANDLES
end


