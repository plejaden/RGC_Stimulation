% how to use: breakpoint in mysimulationbuilder

close all;

indexOffset=size(measurementPoints,2)-3;


xValuesPlot =linspace(1,compartmentNumbers,compartmentNumbers);
xValues=measurementPoints(:,1+indexOffset);
yValues=measurementPoints(:,2+indexOffset);
zValues=measurementPoints(:,3++indexOffset);

plot(xValuesPlot,xValues,xValuesPlot,yValues,xValuesPlot,zValues);
legend('x-values', 'y-values', 'z-values');


