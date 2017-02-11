

if (exist('T','var') && exist('Y','var'))
    
    directoryName = ['output/' 'gatingPlot/' filename(8:end-4)];
    mkdir(directoryName);
    %plotCompartment=53;
    plotCompartment=6;
    
    close all;
    
    hold on;
    plot(T,Y(:,2,plotCompartment), T,Y(:,3,plotCompartment), T, Y(:,4,plotCompartment), T, Y(:,5,plotCompartment));
    legend('h(t)', 'm(t)', 'n(t)', 'c(t)','location','west');

    xlabel('t [ms]');
    ylabel('Magnitude');
    hold off;
    
    
    cleanfigure;
    matlab2tikz([ directoryName '/' 'gating.tikz' ], 'showInfo', false, ...
        'parseStrings',false, ...
        'standalone', false, ...
        'height', '\figureheight', ...
        'width','\figurewidth');
    
else
    display('please load simulation data from output directory');
end