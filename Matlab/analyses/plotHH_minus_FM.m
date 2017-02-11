%% Problem: the resampling is changing the shape of the pulse; its not possible to subtract the pulses without resampling, and with resampling the result is worthless.

%% manual preperation
% step 1: load HH and rename "currentObserver" to "currentObserver_HH"
% step 2: delete other variables from workspace
% step 3: load Fm and rename "currentObserver" to "currentObserver_FM"
% step 4: delete other variables

%% find maximum samples
numberOfSamples=2720;
if size(currentObserver_HH.INa,1)<numberOfSamples
    numberOfSamples=size(currentObserver_HH.INa,1);
end

if size(currentObserver_FM.INa,1)<numberOfSamples
    numberOfSamples=size(currentObserver_FM.INa,1);
end

display(['using ' num2str(numberOfSamples) ' samples']);


%% resample currents

fcn_INa_resampled=resample(currentObserver_FM.INa,numberOfSamples,size(currentObserver_FM.INa,1));
fcn_IK_resampled=resample(currentObserver_FM.IK,numberOfSamples,size(currentObserver_FM.IK,1));
fcn_ICa_resampled=resample(currentObserver_FM.ICa,numberOfSamples,size(currentObserver_FM.ICa,1));
fcn_IL_resampled=resample(currentObserver_FM.IL,numberOfSamples,size(currentObserver_FM.IL,1));
fcn_IKa_resampled=zeros(numberOfSamples,1);

hh_time=currentObserver_HH.time;
time_resampled=linspace(hh_time(1), hh_time(end), numberOfSamples);

hh_INa_resampled=resample(currentObserver_HH.INa,numberOfSamples,size(currentObserver_HH.INa,1));
hh_IK_resampled=resample(currentObserver_HH.IK,numberOfSamples,size(currentObserver_HH.IK,1));
hh_IL_resampled=resample(currentObserver_HH.IL,numberOfSamples,size(currentObserver_HH.IL,1));
hh_ICa_resampled=zeros(numberOfSamples,1);
hh_IKa_resampled=zeros(numberOfSamples,1);

%% getting difference values
iNaDiff=hh_INa_resampled-hh_INa_resampled;
iKDiff=hh_IK_resampled-fcn_IK_resampled;
iKaDiff=hh_IKa_resampled-fcn_ICa_resampled;
iCaDiff=hh_ICa_resampled-fcn_IKa_resampled;
iLDiff=hh_IL_resampled-fcn_IL_resampled;



%% reset values for testing
iNaDiff=zeros(numberOfSamples,1);%fcn_INa_resampled;
iKDiff=hh_INa_resampled;
iKaDiff=zeros(numberOfSamples,1);%fcn_IK_resampled;
iCaDiff=zeros(numberOfSamples,1);%hh_IK_resampled;
iLDiff=zeros(numberOfSamples,1);%hh_INa_resampled + hh_IK_resampled;%

%% ploting
%plot(time_resampled,hh_INa_resampled);
close all;
hold on
plot(hh_time,currentObserver_HH.INa);
plot(time_resampled, iNaDiff, time_resampled,iKDiff, time_resampled, iKaDiff, time_resampled,  iCaDiff, time_resampled,iLDiff);
hold off;
%% export tikz