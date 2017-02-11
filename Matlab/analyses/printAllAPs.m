% load('/home/plejaden/Dokumente/Uni/MasterThesis/trunk/misc/AP_vs_c2c/fcn2010_cat_20deg+-0.04+saiful_AxonIsSoma_epi_SOCD_100cmp_ALLum.mat')
%load('/home/plejaden/Dokumente/Uni/MasterThesis/trunk/misc/AP_vs_c2c/fcn2010_cat_20deg+-0.04+saiful_AxonIsSoma_epi_axon_100cmp_ALLum.mat')

directoryName = 'output/AxonIsSoma_epi_100cmp_ALLum';
mkdir(directoryName);
plotCompartment=53;
%plotCompartment=6;

close all;

hold on;
plot(T_075, Y_075(:,1,plotCompartment),'b-');
plot(T_100, Y_100(:,1,plotCompartment),'b-.');
plot(T_150, Y_150(:,1,plotCompartment),'g-.');
%if exist('Y_150_40','var') plot(T_150_40, Y_150_40(:,1,plotCompartment),'m-'); end
plot(T_200, Y_200(:,1,plotCompartment),'c-.');
plot(T_250, Y_250(:,1,plotCompartment),'m-.');
plot(T_300, Y_300(:,1,plotCompartment),'y-.');
plot(T_400, Y_400(:,1,plotCompartment),'r--');
plot(T_500, Y_500(:,1,plotCompartment),'r-');
%if exist('Y_150_40','var')
%    legend('3.3 um','17.73 um','46.6 um','46.6 um @ -400mV','75.47 um','104.33 um','133.25 um','190.94 um','248.6 um','location','east');
%    axis([0 3.5 -100 60]);
%else
    legend('3.3 um','17.73 um','46.6 um','75.47 um','104.33 um','133.25 um','190.94 um','248.6 um','location','west');
    axis([0 5 -90 60]);
%end


xlabel('t [ms]');
ylabel('Membrane Potential [mV]');
hold off;


cleanfigure;
matlab2tikz([ directoryName '/' 'AP_vs_c2c.tikz' ], 'showInfo', false, ...
    'parseStrings',false, ...
    'standalone', false, ...
    'height', '\figureheight', ...
    'width','\figurewidth');