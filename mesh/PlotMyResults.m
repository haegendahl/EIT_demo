function  PlotMyResults(VisElements, VisNodes,IntpMatrices, sigma,range,figno)

N = numel(VisElements);


clf
minval = range(1);
maxval = range(2);

for ii=1:N
    
    %     minc=min([real(adm(:,ii));imag(adm(:,ii))]);
    %     maxc=max([real(adm(:,ii));imag(adm(:,ii))]);
    %
    figure(ii)
    clf
    axis square
    axis off
    
    H = VisElements{ii};
    g = VisNodes{ii};
    sig = IntpMatrices{ii}*sigma;
    fp = patch('faces',H,'vertices',g,'facevertexcdata',sig,'facecolor','interp','edgecolor','none');
    caxis([minval, maxval])
    
    %set(get(fp,'parent'),'CLim',[minval maxval])
    

    colormap(fireprint())
%    title(['Re[\gamma] at ', 10, num2str(FREQUENCIES(ii),'%.3e'), ' Hz'])
    
    
    
end