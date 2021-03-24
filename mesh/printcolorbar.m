function printcolorbar(cmin, cmax, PRINTFILE)

    fig5=figure(5);
    clf
    set(fig5,'visible','off');
    axis off
    colormap(fireprint())
    colorbar('location','southoutside')
    caxis([cmin cmax]);
    set(fig5,'PaperUnits','centimeters');
    set(fig5,'PaperSize',[20.4800    2.0000]);
    set(fig5,'PaperPosition',[-2.6300   -1.7000   25.0000   16.8100]);
    print(PRINTFILE, '-dpdf');
    %saveas(fig5, printfile, printtype);