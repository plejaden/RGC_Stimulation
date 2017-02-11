clc
T=currentObserver.time;

directoryName = filename(1:end-4);
mkdir(directoryName);


%%
% http://www.doc.ic.ac.uk/~man210/tips/2014/01/22/creating-publication-quality-figures/
% aproach

close all
set(0,'defaulttextinterpreter','none'); % Prevents Matlab from trying to
% interpret latex string.

%% print single currents vs time

iNaAbsolut=0;
iKAbsolut=0;
iKCaAbsolut=0;
iCaAbsolut=0;
iLAbsolut=0;
iStim=0;
T=currentObserver.time;
directoryName = filename(1:end-4);
h.figure = figure;
clf
plotaxis.Xmin=0; plotaxis.Xmax=5; plotaxis.Ymin=-1; plotaxis.Ymax=1;

if  sum(abs(currentObserver.IKA))~=0
    display('fm1997: ion currents vs. time');
    unit = ' mC/cm2';
    
    display('charge per cm2');
    iNaAbsolut=currentObserver.INa; cNa=trapz(T,iNaAbsolut); % [ uA/cm2 * ms = nC/cm2]
    iKAbsolut=currentObserver.IK; cK=trapz(T,iKAbsolut); 
    iKA_Absolut=currentObserver.IKA; cKA=trapz(T,iKA_Absolut); 
    iKCaAbsolut=currentObserver.IKCa; cKCa=trapz(T,iKCaAbsolut);
    iCaAbsolut=currentObserver.ICa; cCa=trapz(T,iCaAbsolut); 
    iLAbsolut=currentObserver.IL; cL=trapz(T,iLAbsolut); 
    iStim=currentObserver.stimulusCurrent; cStim=trapz(T,iStim); 
    
    display(['charge Na: ' num2str(cNa) unit]); display(['charge K: ' num2str(cK) unit ]); display(['charge KA: ' num2str(cKA) unit ]);
    display(['charge KCa: ' num2str(cKCa) unit]);   display(['charge Ca: ' num2str(cCa) unit]); display(['charge L: ' num2str(cL) unit ]);
    
    plot(T, currentObserver.INa /1000, T,currentObserver.IK /1000,         T,currentObserver.IKA / 1000, T, currentObserver.IKCa / 1000, T,  currentObserver.ICa / 1000, T, currentObserver.IL / 1000); % T, currentObserver.stimulusCurrent / 1000
    legend(                'INa',          'IK',                  'IKA',            'IKCa',             'ICa',     'IL','location','east');
    axis([plotaxis.Xmin plotaxis.Xmax plotaxis.Ymin plotaxis.Ymax]);
    xlabel('t [ms]');
    ylabel('i [mA/cm2]');
    
elseif sum(abs(currentObserver.IKCa))~=0
    display('fcn2010: ion currents per cm^2 vs. time');
    unit = ' nC/cm2';
    
    display('charge per cm2');
    iNaAbsolut=currentObserver.INa; cNa=trapz(T,iNaAbsolut);% [ uA/cm2 * ms = nC/cm2]
    iKAbsolut=currentObserver.IK; cK=trapz(T,iKAbsolut);
    % iKA_Absolut=currentObserver.IKA; cKA=trapz(T,iKA_Absolut);
    iKCaAbsolut=currentObserver.IKCa; cKCa=trapz(T,iKCaAbsolut);
    iCaAbsolut=currentObserver.ICa; cCa=trapz(T,iCaAbsolut);
    iLAbsolut=currentObserver.IL; cL=trapz(T,iLAbsolut);
    
    display(['charge Na: ' num2str(cNa) unit]); display(['charge K: ' num2str(cK) unit ]); %display(['charge KA: ' num2str(cKA) ' uC/cm2' ]);
    display(['charge KCa: ' num2str(cKCa) unit]);   display(['charge Ca: ' num2str(cCa) unit]); display(['charge L: ' num2str(cL) unit ]);
    
    plot(T, currentObserver.INa /1000, T,currentObserver.IK / 1000, T, currentObserver.IKCa / 1000, T,  currentObserver.ICa / 1000, T, currentObserver.IL / 1000);
    legend('INa', 'IK', 'IKCa', 'ICa', 'IL','location','east');
    axis([plotaxis.Xmin plotaxis.Xmax plotaxis.Ymin plotaxis.Ymax]);
    xlabel('t [ms]');
    ylabel('i [mA/cm2]');
else
    display('hh: ion currents per cm^2 vs. time');
    unit = ' nC/cm2';
    
    display('charge per cm2');
    iNaAbsolut=currentObserver.INa; cNa=trapz(T,iNaAbsolut); % [ uA/cm2 * ms = nC/cm2]
    iKAbsolut=currentObserver.IK; cK=trapz(T,iKAbsolut);
    iLAbsolut=currentObserver.IL; cL=trapz(T,iLAbsolut);
    
    display(['charge Na: ' num2str(cNa) unit]); display(['charge K: ' num2str(cK) unit]);  display(['charge L: ' num2str(cL) unit]);
    
    plot(T, currentObserver.INa /1000, T,currentObserver.IK / 1000, T, currentObserver.IL / 1000);
    axis([plotaxis.Xmin plotaxis.Xmax plotaxis.Ymin plotaxis.Ymax]);
    lngd=legend('INa', 'IK', 'IL', 'Location', 'east');
    xl=xlabel('t [ms]');
    yl=ylabel('i [mA/cm2]');
end


% save plotted images
cleanfigure;
matlab2tikz([ directoryName '/' 'allCurrentsPerArea.tikz' ], 'showInfo', false, ...
    'parseStrings',false, ...
    'standalone', false, ...
    'height', '\figureheight', ...
    'width','\figurewidth');

%close all;

%% plot absolut currents
T=currentObserver.time;
directoryName = filename(1:end-4);
h1=figure;
membranArea=geometryProperties.compartmentRadius(currentObserver.compartment)*2*pi*geometryProperties.distanceBetrweenCompartments+geometryProperties.compartmentRadius(currentObserver.compartment)^2*pi*2; % [m2]
display(['membranArea: um^2: ', num2str(membranArea * (10^12))]);

iNaAbsolut=0;
iKAbsolut=0;
iKCaAbsolut=0;
iCaAbsolut=0;
iLAbsolut=0;

if  sum(abs(currentObserver.IKA))~=0
    display(' ');
    display('fm1997: ion currents vs. time');
    
    iNaAbsolut=currentObserver.INa * membranArea * 100^2 * 1000;   % [uA /cm2 * cm2 = uA = nA]
    cNa=trapz(T,iNaAbsolut); % [nA * ms = pC ]
    iKAbsolut=currentObserver.IK * membranArea * 100^2 * 1000; cK=trapz(T,iKAbsolut);
    iKA_Absolut=currentObserver.IKA * membranArea * 100^2 * 1000; cKA=trapz(T,iKA_Absolut);
    iKCaAbsolut=currentObserver.IKCa * membranArea * 100^2 * 1000; cKCa=trapz(T,iKCaAbsolut);
    iCaAbsolut=currentObserver.ICa * membranArea * 100^2 * 1000; cCa=trapz(T,iCaAbsolut);
    iLAbsolut=currentObserver.IL * membranArea * 100^2 * 1000; cL=trapz(T,iLAbsolut);
    
    display(['charge Na: ' num2str(cNa) ' pC']); display(['charge K: ' num2str(cK) ' pC' ]); display(['charge KA: ' num2str(cKA) ' pC' ]);
    display(['charge KCa: ' num2str(cKCa) ' pC']);   display(['charge Ca: ' num2str(cCa) ' pC']); display(['charge L: ' num2str(cL) ' pC' ]);
    
    
    plot(T, iNaAbsolut, T,iKAbsolut, T, iKA_Absolut, T, iKCaAbsolut, T,  iCaAbsolut, T, iLAbsolut);
    legend('INa [nA]', 'IK [nA]', 'IKA  [nA]', 'IKCa [nA]', 'ICa [nA]', 'IL [nA]','location','southeast');
    
elseif sum(abs(currentObserver.IKCa))~=0
    display(' ');
    display('fcn2010: ion currents vs. time (absolute)');
    iNaAbsolut=currentObserver.INa * membranArea * 100^2 * 1000; % [uA /cm2 * cm2 = uA = nA]
    cNa=trapz(T,iNaAbsolut); % [nA * ms = pC ]
    iKAbsolut=currentObserver.IK * membranArea * 100^2 * 1000;  cK=trapz(T,iKAbsolut);
    iKCaAbsolut=currentObserver.IKCa * membranArea * 100^2 * 1000;  cKCa=trapz(T,iKCaAbsolut)  ;
    iCaAbsolut=currentObserver.ICa * membranArea * 100^2 * 1000;  cCa=trapz(T,iCaAbsolut)  ;
    iLAbsolut=currentObserver.IL * membranArea * 100^2 * 1000;  cL=trapz(T,iLAbsolut) ;
    
    display(['charge Na: ' num2str(cNa) ' pC']); display(['charge K: ' num2str(cK) ' pC' ]); display(['charge KCa: ' num2str(cKCa) ' pC']);
    display(['charge Ca: ' num2str(cCa) ' pC']); display(['charge L: ' num2str(cL) ' pC' ]);
    
    plot(T, iNaAbsolut , T, iKAbsolut, T, iKCaAbsolut, T, iCaAbsolut, T, iLAbsolut);
    legend('INa [nA]', 'IK [nA]', 'IKCa [nA]', 'ICa [nA]', 'IL [nA]','location','southeast');
else
    display(' ');
    display('hh: ion currents vs. time (absolute)');
    iNaAbsolut=currentObserver.INa * membranArea * 100^2 * 1000; % [uA /cm2 * cm2 = uA = nA]
    cNa=trapz(T,iNaAbsolut); % [nA * ms = pC ]
    iKAbsolut=currentObserver.IK * membranArea * 100^2 * 1000;  cK=trapz(T,iKAbsolut); 
    iLAbsolut=currentObserver.IL * membranArea * 100^2 * 1000; cL=trapz(T,iLAbsolut) ;
    
    display(['charge Na: ' num2str(cNa) ' pC']); display(['charge K: ' num2str(cK) ' pC' ]); display(['charge L: ' num2str(cL) ' pC' ]);
    
    plot(T, iNaAbsolut , T, iKAbsolut, T, iLAbsolut);
    legend('INa [nA]', 'IK [nA]', 'IL [nA]','location','southeast');
    xlabel('t [ms]');
    ylabel('Q [pC]');
end



% save plotted images
cleanfigure;
matlab2tikz([ directoryName '/' 'allCurrentsAbsolute.tikz' ], 'showInfo', false, ...
    'parseStrings',false, ...
    'standalone', false, ...
    'height', '\figureheight', ...
    'width','\figurewidth');



%close all;


% %% plot distributuin of currents
% h2=figure;
% membranArea=geometryProperties.compartmentRadius(currentObserver.compartment)*100*2*pi*geometryProperties.distanceBetrweenCompartments*100;
%
% if  sum(abs(currentObserver.IKA))~=0
%     display('fm1997: electrical charge of the ion channels');
%     iNaAbsolut=currentObserver.INa * membranArea * 10000; cNa=trapz(T,iNaAbsolut); % [ ms * uA/cm2 * ms = uC]
%     iKAbsolut=currentObserver.IK * membranArea * 10000; cK=trapz(T,iKAbsolut); % [ ms * nA = uC]
%     iKA_Absolut=currentObserver.IKA * membranArea * 10000; cKA=trapz(T,iKA_Absolut); % [ ms * nA = uC]
%     iKCaAbsolut=currentObserver.IKCa * membranArea * 10000; cKCa=trapz(T,iKCaAbsolut); % [ ms * nA = uC]
%     iCaAbsolut=currentObserver.ICa * membranArea * 10000; cCa=trapz(T,iCaAbsolut); % [ ms * nA = uC]
%     iLAbsolut=currentObserver.IL * membranArea * 10000; cL=trapz(T,iLAbsolut); % [ ms * nA = uC]
%     charge=[ abs(cNa) abs(cK) abs(cKA) abs(cKCa) abs(cCa) abs(cL) ]*1000; % [ uC * 1000 = pC]
%
%     display(['charge Na: ' num2str(cNa) ' pC']); display(['charge K: ' num2str(cK) ' pC' ]); display(['charge KA: ' num2str(cKA) ' pC' ]);
%     display(['charge KCa: ' num2str(cKCa) ' pC']);   display(['charge Ca: ' num2str(cCa) ' pC']); display(['charge L: ' num2str(cL) ' pC' ]);
%
%     labels = {'INa','IK','IKA','IKCa','ICa','IL'};
%     pie(charge,labels)
%
% elseif sum(abs(currentObserver.IKCa))~=0
%     display('fcn2010: electrical charge of the ion channels');
%     iNaAbsolut=currentObserver.INa * membranArea * 1000; cNa=trapz(T,iNaAbsolut); % [ ms * nA = pC]
%     iKAbsolut=currentObserver.IK * membranArea * 1000; cK=trapz(T,iKAbsolut); % [ ms * nA = pC]
%     iKCaAbsolut=currentObserver.IKCa * membranArea * 1000; cKCa=trapz(T,iKCaAbsolut); % [ ms * nA = pC]
%     iCaAbsolut=currentObserver.ICa * membranArea * 1000; cCa=trapz(T,iCaAbsolut); % [ ms * nA = pC]
%     iLAbsolut=currentObserver.IL * membranArea * 1000; cL=trapz(T,iLAbsolut); % [ ms * nA = pC]
%     charge=[ abs(cNa) abs(cK) abs(cKCa) abs(cCa) abs(cL) ]*1000; % [ ms * nA = pC]
%

%
%     labels = {'INa','IK','IKCa','ICa','IL'};
%     pie(charge,labels)
% else
%
%     display('hh: electrical charge of the ion channels');
%     iNaAbsolut=currentObserver.INa * membranArea * 1000; cNa=trapz(T,iNaAbsolut); % [ ms * nA = pC]
%     iKAbsolut=currentObserver.IK * membranArea * 1000; cK=trapz(T,iKAbsolut); % [ ms * nA = pC]
%     iLAbsolut=currentObserver.IL * membranArea * 1000; cL=trapz(T,iLAbsolut); % [ ms * nA = pC]
%     charge=[ abs(cNa) abs(cK) abs(cL) ]*1000; % [ ms * nA = pC]
%
%     display(['charge Na: ' num2str(cNa) ' pC']); display(['charge K: ' num2str(cK) ' pC' ]); display(['charge L: ' num2str(cL) ' pC' ]);
%
%     labels = {'INa','IK','IL'};
%     pie(charge,labels)
%
% end
%
% % save plotted images
% matlab2tikz([ directoryName '/' 'pieChartCurrents.tikz' ], 'showInfo', false, ...
%     'parseStrings',false, ...
%     'standalone', false, ...
%     'height', '\figureheight', ...
%     'width','\figurewidth');
%
% close all;
% clear all;
