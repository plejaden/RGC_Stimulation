function showStimMap( geometryDistance , plotElectrodes)

close all;

resolution=10; % um
electrodeDiameter = 40; %um
centerX=1.5; % mm %center (y axis) of the center electrode
centerY=2; % mm %center (y axis) of the center electrode


cd /home/plejaden/Dokumente/Uni/MasterThesis/trunk;

if(geometryDistance==075)
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_075um_SOCB_res10_steps0.1_t50_mr300_mp2500_solutions.csv');
elseif(geometryDistance==200)
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_200um_SOCB_res10_steps0.1_t50_mr300_mp2500_solutions.csv');
elseif(geometryDistance==215)
    %200 um geometry with ode15
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_223um_SOCB_res10_steps0.1_t25_mr60_mp2500_solutions.csv');
elseif(geometryDistance==223)
    %200 um geometry with ode23
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_223um_SOCB_res10_steps0.1_t25_mr60_mp2500_solutions.csv');
elseif(geometryDistance==224)
    %200 um geometry with ode23s
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_224um_SOCB_res10_steps0.1_t25_mr60_mp2500_solutions.csv');
elseif(geometryDistance==245)
    %200 um geometry with ode45
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_245um_SOCB_res10_steps0.1_t25_mr60_mp2500_solutions.csv');
elseif(geometryDistance==500)
    stimValues=csvread('./misc/stimulationMap_raw/10um/saiful_SomaIsAxon_epiretinal_500um_SOCB_res10_steps0.1_t50_mr300_mp2500_solutions.csv');
else
    display('ERROR: DATA FOR GEOMETRY NOT AVAILABLE');
    return;
end

%value 0 if error
%clearedData=(stimValues(:,3).*(stimValues(:,4)-1)).*-1;
x=stimValues(:,1);
y=stimValues(:,2);
dataPoints=stimValues(:,3);

minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);

sizex=size(unique(stimValues(:,1)),1);
sizey=size(unique(stimValues(:,2)),1);

%% for development only: mark specific point and print coords
% testingPoint=1;
% dataPoints(testingPoint)=2;
% display(['x-value:',num2str(stimValues(testingPoint,1))]);
% display(['y-value:',num2str(stimValues(testingPoint,2))]);

%% for development only: mark first point
%dataPoints(1)=2; stimValues(1,4)=0;

%% for development only: mark specific coords and print point
% needleX=maxx;
% needleY=0.0019725;
% i=1;
% for i= 1:size(stimValues,1)
%     if stimValues(i,1)==needleX && stimValues(i,2)==needleY
%         dataPoints(i)=2;
%         display(['i-value:',num2str(i)]);
%         display(['x-value:',num2str(stimValues(i,1))]);
%         display(['y-value:',num2str(stimValues(i,2))]);
%
%     end % if stimValues(1,i)==needlex && stimValues(2,i)==needleY
% end % for i= 1:size(stimValues,1)


% find y position
needleX=minx;
needleY=centerY/1000;
scaleYpos=1;
scaleYposFound=false;
for i= 1:size(stimValues,1)
    if stimValues(i,1)==needleX && stimValues(i,2)>=needleY-(resolution/1000/1000)/2 && ~scaleYposFound
        scaleYposFound=true;
        display(['Y-center-position:',num2str(i)]);
        display(['x-value:',num2str(stimValues(i,1))]);
        display(['y-value:',num2str(stimValues(i,2))]);
        scaleYpos=i;
%        dataPoints(i)=3; stimValues(i,4)=0;
    end % if stimValues(1,i)==needlex && stimValues(2,i)==needleY
end % for i= 1:size(stimValues,1)


heatMapImageLeft = zeros(sizey,sizex);

k=0;
for thisY=sizey:-1:1
    for thisX=1:sizex
        k=k+1;
        
        thisDatapoint = dataPoints(k);
        if (stimValues(k,4)>0)
%            thisDatapoint=NaN;
        end
        
        heatMapImageLeft(thisY,thisX) = thisDatapoint;
        
    end % for thisY=1 : sizey
    
end % for thisX=1 : sizex


heatMapImageRight=fliplr(heatMapImageLeft);
%heatMapImageRight=zeros(size(heatMapImageLeft)); % disable for debugging
heatMapImage = [heatMapImageLeft,heatMapImageRight];

colormap('jet');
imagesc(heatMapImage);
xlabel('Pixel: x axis') % x-axis label
ylabel('Pixel: y axis') % y-axis label
axis image;
c=colorbar;
set(c, 'ylim', [-5 0])
ylabel(c,'Activation potential in V') ;

if( plotElectrodes)
    %% add shape of electrodes
    r=(electrodeDiameter/2)/resolution;
    
    scaleY=sizey-(((scaleYpos-1)/sizex))+0.5; %perfect centered
    scaleX=(maxx-minx)/(resolution/1000/1000) +1 +1/2; %perfect centered
    
    normalDistA = ((geometryDistance/resolution)/2);
    normalDistC=normalDistA/sin(60*2*pi/360);
    normalDistB=normalDistC*0.5;
    
    hold on;
    % center
    plotCircle( r, scaleX,scaleY ); % center
    plotCircle( r, scaleX+(-normalDistB), scaleY+(normalDistA) ); % Guard I
    plotCircle( r, scaleX+(+normalDistB), scaleY+(normalDistA) ); % Guard II
    plotCircle( r, scaleX+(normalDistC), scaleY+(0) ); % Guard III
    plotCircle( r, scaleX+(+normalDistB), scaleY+(-normalDistA) ); % Guard IV
    plotCircle( r, scaleX+(-normalDistB), scaleY+(-normalDistA) ); % Guard V
    plotCircle( r, scaleX+(-normalDistC), scaleY+(0) ); % Guard VI
    
    hold off
end % if( plotElectrodes)
end

