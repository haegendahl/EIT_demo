function  PrintMyResults(VisElements, VisNodes,IntpMatrices, sigma, slices,range, pathtarget)

if ~exist(pathtarget,'dir')
    mkdir(pathtarget);
end

printtype = 'png';

N = numel(VisElements);


minval = range(1);
maxval = range(2);

for ii=1:N
    
    fig_number = ['f', num2str(ii)];
    f = figure(ii);

    clf
    axis square
    axis off
    H = VisElements{ii};
    g = VisNodes{ii};
    sig = IntpMatrices{ii}*sigma;
    fp = patch('faces',H,'vertices',g,'facevertexcdata',sig,'facecolor','interp','edgecolor','none');
    caxis([minval, maxval])
 
    printfile = [pathtarget '/' slices '_' fig_number '.png'];
    saveas(f, printfile, printtype);
    close(f)
    
    
    %set(get(fp,'parent'),'CLim',[minval maxval])
    
    
    %colormap(fireprint())
    %    title(['Re[\gamma] at ', 10, num2str(FREQUENCIES(ii),'%.3e'), ' Hz'])
    
    
    
end