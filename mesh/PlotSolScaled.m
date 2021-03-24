function [fh,fp]=PlotSolScaled(g,H,V,minval,maxval);

% Modified by M. Vauhkonen 19.9.1994 from version paraevol.m.
% Modified by M. Aku Seppänen 22.11.1999 from version plotsol.m.

[VN,VM]=size(V);
[gN,gM]=size(g);
[HN,HM]=size(H);

for ii=1:VM
  view(2);
  hold on
%  set(fh(ii),'defaultpatchfacecolor',[0 0 0],'defaultpatchedgecolor',[1 1 1]);
  for ij=1:HN
    X=g(H(ij,:),1);
    Y=g(H(ij,:),2);
    Z=V(H(ij,:),ii);
    %%%%%%%%%%fp=patch(X,Y,Z,Z);
    fp=patch(X,Y,Z);
%     fp=patch(X,Y,Z,'b');
%    set(fp,'zdata',Z,'facecolor',[0 0 0],'edgecolor',[1 1 1]);
    set(fp,'edgecolor','none');
    hold on
  end
  set(get(fp,'parent'),'CLim',[minval maxval])  
end
axis('square')
%shading interp
                                       
