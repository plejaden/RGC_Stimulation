function [] = sliderPlot(T,Y,geometryProperties)
% Plot different plots according to slider location.

S.T=T;
S.Y=Y;
S.comNumber=size(Y,3); % number of compartments to plot
S.retinaModel=geometryProperties.model;  % the comsol geometry for retina
S.comNumber=geometryProperties.numberOfCompartments;
S.ev=geometryProperties.externalVoltages;
S.compNeuronType=geometryProperties.compFeatures;
S.filename=geometryProperties.filename;


% add aditional artificial compartments for:
S.comNumber=S.comNumber+1; % overview of all voltages (slide 2)
S.comNumber=S.comNumber+1; % voltage over copartment
S.comNumber=S.comNumber+1; % current sum

% S.fh = figure('units','pixels','position',[100 100 800 500],'name','slider_plot','numbertitle','off','resize','off','toolbar','figure');
S.fh = figure('units','pixels','position',[200 200 1600 1000],'name','slider_plot','numbertitle','off','resize','off','toolbar','figure');

plotFirstWindow(S);

silderSteps=1/(S.comNumber);

S.sl = uicontrol('style','slide','unit','pix','position',[20 10 1520 30],'min',1,'max',(S.comNumber),'val',1,'sliderstep',[silderSteps silderSteps],'callback',{@sl_call,S});
%S.sl = uicontrol('style','slide','unit','pix','position',[20 10 760 30],'min',1,'max',(S.comNumber),'val',1,'sliderstep',[silderSteps silderSteps],'callback',{@sl_call,S});
end %function


function [] = sl_call(varargin)
%% 100 compartments
sliderPlotMultiplicator = 1.0049;
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
sliderPosition=get(h,'value');

actualCompartment=round(sliderPosition*sliderPlotMultiplicator)-1;

if (actualCompartment <0)
    actualCompartment=0;
end

% first slider
if (actualCompartment==0)
    plotFirstWindow(S);
    return;
end

% second slider
if (actualCompartment==1)
    plotSecondWindow(S,actualCompartment);
    return;
end

%pre-last slider
if (actualCompartment==(S.comNumber-1))
    plotExternalPotential(S);
    return;
end

% last slider
if (actualCompartment>=(S.comNumber))
    plotAbsCurrentSum(S);
    return;
end

% plot the real compartments,
% plotPotential(str2num(actualCompartment)-1,S);
% plotGates(str2num(actualCompartment)-1,S);
plotPotential(actualCompartment-1,S);
plotGates(actualCompartment-1,S);

end %function

function []=plotPotential(compartment,S)
T=S.T; Y=S.Y; %retinaModel=S.retinaModel;
if  (isempty (S.compNeuronType))
    display(['showing compartment: ' num2str(compartment) ]);
else
    display(['showing compartment: ' num2str(compartment) ' (' char(S.compNeuronType{compartment}) ')']);
%    display(['showing compartment: ' num2str(compartment) ' (' cell2mat(S.compNeuronType{compartment}) ')']);
end
%find automatic and static axis scaling for all compartments
axisX.min=min(T);
axisX.max=max(T);
axisY.max=max(max(Y(:,1,:)))*1.1;
if axisY.max < 0
    axisY.max=0;
end
axisY.min=min(min(Y(:,1,1)))*1.2;
%compartment=compartment+1;

subplot(2,1,1);
plot(T,Y(:,1,compartment));
axis([axisX.min axisX.max axisY.min axisY.max]);
xlabel('t [ms]'); ylabel('membran potential [mV]');
if  (isempty (S.compNeuronType))
    title(['Membrane potential of compartment ' num2str(compartment) ' vs time']);
else
    title(['Membrane potential of compartment ' num2str(compartment) ' vs time' ' (' char(S.compNeuronType{compartment}) ')']);
end

end


function []=plotGates(compartment,S)
T=S.T; Y=S.Y; %retinaModel=S.retinaModel;
subplot(2,1,2);
xlabel('t [ms]'); %ylabel('V(t)');
title(['Gating variables of compartment ' num2str(compartment) ' vs time']);
if (size(Y,2) == 4)
    %disp('no p channel')
    plot(T,Y(:,2), T,Y(:,3), T, Y(:,4));
    legend('h(t)', 'm(t)', 'n(t)');
elseif (size(Y,2) == 5)
    %disp('p channel')
    %    plot(T,Y(:,2), T,Y(:,3), T, Y(:,4), T, Y(:,5));
    plot (T,Y(:,2,compartment), T,Y(:,3,compartment), T, Y(:,4,compartment), T, Y(:,5,compartment));
    legend('h(t)', 'm(t)', 'n(t)', 'p(t)');
elseif (size(Y,2) == 6)
    %disp('p channel')
    plot(T,Y(:,2,compartment), T,Y(:,3,compartment), T, Y(:,4,compartment), T, Y(:,5,compartment));
    legend('h(t)', 'm(t)', 'n(t)', 'c(t)');
elseif (size(Y,2) == 8)
    %disp('p channel')
    plot(T,Y(:,2,compartment), T,Y(:,3,compartment), T, Y(:,4,compartment), T, Y(:,5,compartment), T, Y(:,6,compartment), T, Y(:,7,compartment));
    legend('h(t)', 'm(t)', 'n(t)', 'c(t)', 'a(t)', 'hA(t)');
end
%title('Gating variables vs time');

end %function


function []=plotFirstWindow(S)
retinaModel=S.retinaModel;
subplot(1,1,1);

if ischar(retinaModel) || isempty(retinaModel)
    display('No geometry information found - Plotting first compartment');
    randomCompartment=1; % TODO: repace with real random fnk, but for developing reasons, I use 1 %
    plotPotential(randomCompartment,S);
    plotGates(randomCompartment,S);
    
else
    display('showing gemometry plot');
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    mphgeom(retinaModel,'geom1','Facealpha',0.7);
end
end %function

function []=plotSecondWindow(S,~)
subplot(1,1,1);
T=S.T; Y=S.Y; %retinaModel=S.retinaModel;
iterations=size(Y);
iterations=iterations(1,3);

%% comment in for propergation diagram
close all;
%iterations=20;
display(['showing compartments: 1-' num2str(iterations) ]);

compartmentoffset=50;
%plot(T,0); % plot zero-line
color = [ repmat({'red'},1,1) ; repmat({'black'},3,1); repmat({'green'},3,1); repmat({'blue'},93,1) ];

hold on;
axis auto;
for i=1:iterations;
    plot(T,Y(:,1,i)+(i-1)*compartmentoffset,color{i});
end;
xlabel('t [ms]');
ylabel(['membran potential [mV], compartments seperated by ',num2str(compartmentoffset) ]);
title('Membrane potential all compartments vs time [ms]');
hold off;

end

function []=plotExternalPotential(S)
subplot(1,1,1);
title('voltage [mV] over compartment');
plot(S.ev,'DisplayName','Voltages in mV','YDataSource','Voltages in [mV]');
end

function []=plotAbsCurrentSum(S)
global activation;

% avoid a crash if this function is started with only stored values
if (isempty(activation))
    % simply print first compartment if no current sum is available
    plotPotential(1,S);
    plotGates(1,S);
    return;
end

subplot(1,1,1);
display('showing absolute current' );
title('energy [mA] over compartment');
p=plot(activation.current./10^6,'DisplayName','activation current','YDataSource','Current in [A]');
set(p, 'Color', 'green');
return;
end