function PlotScaledSlices(slices, pathtarget,type_of_slice,print)

close all

printtype = 'png';

for i= 1 : numel(slices)
    clf
    fig_number = ['f', num2str(i)];
    f = figure(i);

    %PlotSolScaled(slices(slice_num ).g,slices(slice_num ).H,slices(slice_num ).theta,0.5,2);   
    axis equal
    axis off
    patch('faces',slices(i).H, 'vertices',slices(i).g, 'facevertexcdata',slices(i).theta,'facecolor','interp','edgecolor','none');
    colormap(fireprint());
    caxis([0.5,2.1])
    
    if(strcmp(print,'yes') ||strcmp(print,'Yes') ||strcmp(print,'YES') )    
    printfile = [pathtarget '/' type_of_slice '_' fig_number '.png'];
    saveas(f, printfile, printtype);
    end
    close(f)
    
end
