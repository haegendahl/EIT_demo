function PlotGradientSlices(slices, pathtarget,type_of_slice,print)

close all

printtype = 'png';

for i= 1 : numel(slices)
    
    fig_number = ['f', num2str(i)];
    f = figure(i);

    [R Ai] = constructTVIntMat(slices(i).g,slices(i).H);
    PlotGradient(slices(i).theta, R,slices(i).g, slices(i).H, i)
    colormap(fireprint());
    if(strcmp(print,'yes') ||strcmp(print,'Yes') ||strcmp(print,'YES') ) 
    printfile = [pathtarget '/' 'grad_' type_of_slice '_' fig_number '.png'];
    saveas(f, printfile, printtype);
    end
    close(f)
    
end