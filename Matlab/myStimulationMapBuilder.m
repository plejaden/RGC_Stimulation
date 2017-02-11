classdef myStimulationMapBuilder
    %% 2015-02-28 11:29:58
    
    %% Start with:
    % clear; myStimulationMapBuilder('fm1997','-','2electrodes_epiretinal');
    methods
        % class constructor
        function this = myStimulationMapBuilder(simulationModelName, searchDirection, electrodeConfiguration)
            compartmentNumbers= 100;
            potentialMultiplicator=0;
            
            mapResolution=10;%in um %10
            potentialSteps=0.1;% in V => 1mV %0.001
            timeout=300;%s defined in mathPlugin files.
            maxRuns=50; % 50
            maxPixel=2500; % 2500
            
            
            ERRORCODE_MAXRUNS=1;
            ERRORCODE_MAXPIXELS=2;
            ERRORCODE_TIMEOUT=3;
            ERRORCODE_EXCEPTION=4;
            
            
            cd /home/plejaden/Dokumente/Uni/MasterThesis/trunk;
            
            solutionFileName=[electrodeConfiguration,'_res',num2str(mapResolution),'_steps',num2str(potentialSteps),'_t',num2str(maxRuns),'_mr',num2str(timeout),'_mp',num2str(maxPixel),'_solutions.csv'];
            
            global error;
            error = ''; % no error occured
            starttime=datestr(now, 'yyyy-mm-dd HH:MM:SS');
            display(' ');
            display([ 'Simulation started at ' starttime]);
            close all;
            warning ('off','all');
            
            tic;
            numberOfParameters=nargin();
            
            if(numberOfParameters ~= 3 )
                display('Quick-guide for using this system');
                display('myStimulationMapBuilder(''fcn2010_cat'',''-'',''2electrodes_epiretinal'')');
                display(' ');
                display('#1 parameter: used model for calculation. Implemented models are:');
                dir('membraneModels/*.m')
                display(' ');
                display('#2 parameter: voltage signum. - for search in neg, + for search in pos potential');
                display(' ');
                display('#3 parameter: Geometry. Implemented geometries are:');
                this.displayPhysicModels();
                return;
            end
            
            display('available Geometries:');
            this.displayPhysicModels();
            
            if(numberOfParameters == 3)
                
                % create retina, nerve and electrode geometry, mesh, run physics %
                if (isempty (error))
                    display('## calculate finite element potentials ##'); toc;
                    geometryProperties=this.organizeModel(this,electrodeConfiguration,compartmentNumbers,potentialMultiplicator);
                    geometryProperties.filename=[ simulationModelName,'+', num2str(potentialMultiplicator),'+', electrodeConfiguration,'+',num2str(compartmentNumbers),'+',datestr(now, 'yyyy-mm-dd_HH:MM:SS') ];
                    geometryProperties.mapResolution=mapResolution;
                    [xMin,yMin] = this.findStartCoords(geometryProperties);
                    [xMax,yMax] = this.findEndCoords(geometryProperties);
                end
                
                allRunCount=0;
                pixel=0;
                geometryProperties = this.adoptNeuronCoords(geometryProperties,xMin,yMin,xMax,yMax,pixel);
                % correction of X.points if needed
                %pixel=92; geometryProperties.xNeuronStart=0.0012313; geometryProperties.yNeuronStart=0.00183;
                while(~this.calculationComplete(geometryProperties,xMin,yMin,pixel) && pixel < maxPixel)
                    display(' '); display(' '); display(' '); display(' '); display(' '); display(' '); display(' '); display(' ');
                    display(['%%##%%##%%##%%##%%##']);
                    display(['# pixel: ',num2str(pixel)]);
                    display(['# corods x: ',num2str(geometryProperties.xNeuronStart) , ' xMax: ',num2str(xMax)]);
                    display(['# corods y: ',num2str(geometryProperties.yNeuronStart) , ' yMax: ',num2str(yMax)]);
                    toc;
                    
                    
                    %% find lowest threshold for x-y point
                    run=0;
                    solutionFoundForSpot = false;
                    while(~solutionFoundForSpot)
                        thistime=datestr(now, 'yyyy-mm-dd HH:MM:SS');
                        display(' ');
                        display(['%%##%%##%%##%%##%%##']);
                        display(['# run: ',num2str(run) ,' (' ,num2str(allRunCount),')']);
                        display(['# point: ',num2str(pixel)]);
                        display(['# potential: ',num2str(potentialMultiplicator),'V']);
                        display(['# time: ' thistime]);
                        toc;
                        
                        geometryProperties = this.getExVoltages(geometryProperties);
                        
                        % envokes the nerve model for calculating the ion channel behaviour
                        if (isempty (error))
                            simulationModel = this.envokeSimulationModel(simulationModelName, geometryProperties);
                        end
                        
                        
                        try
                            %Calculate the model
                            if (isempty (error))
                                display(' '); display('## solveing ode ##'); toc;
                                [T, Y] = simulationModel.run(geometryProperties);
                            end
                        catch
                            
                            display('Error: No solution detected');
                            dlmwrite(solutionFileName,[geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator,ERRORCODE_EXCEPTION],'delimiter',',','-append');
                            solutionFoundForSpot=true;
                            Y=zeros(10);
                            T=zeros(10);
                        end
                        
                        
                        if (isempty (error))
                            
                            % check if there are some differences in the resulting matrix, which may be an AP %
                            if (~this.findAPinData(Y))
                                display(['No AP spike detected! (run ',num2str(run),')']); toc;
                                
%                                 if(mod(run,25)==0)
%                                     this.plotCheck(Y,T,geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator);
%                                 end
                                
                                if(strfind(searchDirection, '-'))
                                    potentialMultiplicator = potentialMultiplicator - potentialSteps;
                                else
                                    potentialMultiplicator = potentialMultiplicator + potentialSteps;
                                end
                                geometryProperties.potMultiplicator=potentialMultiplicator;
                                
                                
                            else
                                global currentObserver;
                                
                                solutionFoundForSpot=true;
                                dlmwrite(solutionFileName,[geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator],'delimiter',',','-append');
                                
                                this.plotCheck(Y,T,geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator);
                                
                            end
                        end
                        
                        if(~solutionFoundForSpot && sum(T(:))==0)
                            solutionFoundForSpot=true;
                            dlmwrite(solutionFileName,[geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator,ERRORCODE_TIMEOUT],'delimiter',',','-append');
                            display('ERRORCODE_TIMEOUT');
                        end
                        
                        
                        if (~solutionFoundForSpot && run>=maxRuns)
                            display('No solution detected');
                            dlmwrite(solutionFileName,[geometryProperties.xNeuronStart,geometryProperties.yNeuronStart,potentialMultiplicator,ERRORCODE_MAXRUNS],'delimiter',',','-append');
                            solutionFoundForSpot=true;
                        end
                        
                        geometryProperties.potMultiplicator=potentialMultiplicator;
                        run=run+1;
                        allRunCount=allRunCount+1;
                    end %  while(~solutionFoundForSpot)
                    
                    pixel=pixel+1;
                    display(['solution for point x:',num2str(geometryProperties.xNeuronStart),' y:',num2str(geometryProperties.yNeuronStart),' found:',num2str(potentialMultiplicator),'V']);
                    potentialMultiplicator=0;
                    geometryProperties = this.adoptNeuronCoords(geometryProperties,xMin,yMin,xMax,yMax,pixel);
                    geometryProperties.potMultiplicator=potentialMultiplicator;
                    
                end % for(i=0; calculationComplete(geometryProperties,xMin,yMin,iteration);i=i+1)
                
                
                
            elseif (numberOfParameters > 4)
                %disp('checkpoint');
                fprintf('%d parameters, which is too much (3 parameters needed: Model, impulse Type and electrodeConfiguration)',nargin);
            end %if (nargin == 3)
            display(' '); toc;
            display(' ');
            display(' ');
        end %function myModelBuilder
        
    end %methods
    
    
    methods (Static)
        function geometryProperties=organizeModel(this,modelType,compartmentNumbers,potentialMultiplicator)
            
            
            geometryProperties.potMultiplicator=potentialMultiplicator;
            geometryProperties.numberOfCompartments=compartmentNumbers;
            geometryProperties.iStim=0;
            geometryProperties.rho = 70.6/1000; % 70.6 Ohm/um => Ohm/m
            geometryProperties.rhoLayer=18.2/1000; % 18.2 Ohm/um => Ohm/m
            
            %% Values, generated without geometry.
            % stimulate with already measured values
            
            if strfind(modelType, 'prestored')
                display(' %% Using stored geometry results');
                load (['retinaElectrodeModels/' 'preRenderedModels/' modelType '.mat'])
                geometryProperties.externalVoltages=geometryProperties.externalVoltages*potentialMultiplicator;
                geometryProperties.compFeatures = [ repmat({'Soma'},1,1) ; repmat({'IS'},3,1); repmat({'SOCB'},3,1); repmat({'Axon'},93,1) ];
                geometryProperties.compartmentRadius=[ repmat(10*10^-6,1,1) ; repmat(10^-6,3,1); repmat(10^-6,3,1); repmat(10^-6,93,1) ];
                return;
            end
            
            
            %stimulate with dirac
            if strcmpi(modelType,'testPointStimulation')
                display('TEST-MODE: Using Point function');
                geometryProperties.model='noModel-precalculated Values only';
                geometryProperties.numberOfCompartments=compartmentNumbers;
                geometryProperties.electrodeDistance=1e-04;
                geometryProperties.compartmentRadius=5e-03; % working, but not correct
                geometryProperties.compartmentRadius=1e-05; % = 10um % usually 0.2-20um
                geometryProperties.nerveLength=3e-3;
                geometryProperties.distanceBetrweenCompartments=geometryProperties.nerveLength/compartmentNumbers; %m
                
                geometryProperties.compartmentRadius=5e-05; %m
                
                oneComp=zeros(geometryProperties.numberOfCompartments,1);
                oneComp(6)=1;
                geometryProperties.externalVoltages=oneComp;
                geometryProperties.externalVoltages=geometryProperties.externalVoltages*1000; %convert volts to mV
                
                return;
            end
            
            % no external stimulation
            if strcmpi(modelType,'noStimulation')
                display('TEST-MODE: Disabled all stimulus');
                geometryProperties.model='noModel-precalculated Values only';
                geometryProperties.rho = 0.0706;
                geometryProperties.numberOfCompartments=compartmentNumbers;
                
                geometryProperties.compartmentRadius=[ repmat(10^-6,100,1) ];
                geometryProperties.electrodeDistance=6*10^-5;
                geometryProperties.nerveLength=3e-3;
                geometryProperties.distanceBetrweenCompartments=geometryProperties.nerveLength/compartmentNumbers; %m
                
                geometryProperties.externalVoltages=zeros(1,geometryProperties.numberOfCompartments);
                return;
                
            end
            
            
            % no external stimulation
            if strcmpi(modelType,'intercellularStimulationOnly')
                display('_-_-_-_-_-_-')
                display('INTERCELLULAR-STIMULATION:');
                display('Disabled all stimulus but the "potentialMultiplicator" [ua/cm^2] at the first compartment');
                geometryProperties.model='noModel-virtual soma stimulaton only';
                geometryProperties.rho = 50.5;
                geometryProperties.numberOfCompartments=compartmentNumbers;
                geometryProperties.iStim=potentialMultiplicator;
                
                geometryProperties.compFeatures = [ repmat({'Soma'},1,1) ; repmat({'IS'},3,1); repmat({'SOCB'},3,1); repmat({'Axon'},93,1) ];
                geometryProperties.compartmentRadius=[ repmat(10*10^-6,1,1) ; repmat(10^-6,3,1); repmat(10^-6,3,1); repmat(10^-6,93,1) ];
                geometryProperties.electrodeDistance=6*10^-5;
                geometryProperties.nerveLength=3e-3;
                geometryProperties.distanceBetrweenCompartments=geometryProperties.nerveLength/compartmentNumbers; %m
                
                geometryProperties.externalVoltages=zeros(1,compartmentNumbers);
                return;
                
            end
            
            %% Usage of real geometry
            geoFuncName=['geometry_' modelType ];
            gdescName=['description_' modelType ];
            
            if (exist(geoFuncName, 'file') && exist(gdescName, 'file') )
                gmodel=str2func(geoFuncName);
                gdesc=str2func(gdescName);
            else
                display('Geometry: no correct geometry specified, using a very simple geometry');
                gmodel=simple2ElectrodeConfiguration();
                gdesc=description_simple2ElectrodeConfiguration();
            end
            
            % envoking the gmodel and gdesc functionpointers.
            geometryProperties.model=gmodel();
            geometryProperties.desc=gdesc(); toc;
            
            geometryProperties.numberOfCompartments=compartmentNumbers;
            geometryProperties=this.readFEMModel(geometryProperties); % add most of the properties
            
        end % geometryProperties=organizeModel(this,modelType,compartmentNumbers,potentialMultiplicator)
        
        %% get the external Voltages form the geometry
        function geometryProperties=readFEMModel(geometryProperties)
            
            retinaModel=geometryProperties.model;
            compartmentNumbers=geometryProperties.numberOfCompartments;
            geometryProperties.distanceBetrweenCompartments=10.69*10^-5; % ??
            
            import com.comsol.model.*
            import com.comsol.model.util.*
            
            %% Step 1: get the names of all retinal geometries
            nerveGeometryNames=cellstr( geometryProperties.desc.neuronPartNames(:,1) );
            electrodeGeometryNames=geometryProperties.desc.electrodeGeometryNames;
            
            % get variable definition, if there are any
            if max(strcmp(fieldnames(geometryProperties.desc),'variables'))
                for i=1:length(geometryProperties.desc.variables)
                    evalc(char(geometryProperties.desc.variables(i)));
                end
            end
            
            
            %% Step 2: get the coordinates from the geometries
            numberOfNeuronparts= size(nerveGeometryNames(:,1),1);
            neuronCoords.x=zeros(numberOfNeuronparts,1);
            neuronCoords.y=zeros(numberOfNeuronparts,1);
            neuronCoords.z=zeros(numberOfNeuronparts,1);
            neuronCoords.h=zeros(numberOfNeuronparts,1);
            neuronCoords.r=zeros(numberOfNeuronparts,1);
            neuronCoords.name=cell(numberOfNeuronparts,1);
            
            for i=1:numberOfNeuronparts
                neuronPartName=nerveGeometryNames(i,:);
                propertiesNerve=mphgetproperties(retinaModel.geom('geom1').feature(neuronPartName));
                display(['reading neuronal part: ',neuronPartName{1}]);
                neuronCoords.name{i}=neuronPartName;
                neuronCoords.x(i)=eval(propertiesNerve.x);
                neuronCoords.y(i)=eval(propertiesNerve.y);
                neuronCoords.z(i)=eval(propertiesNerve.z);
                %                  display(['__x:' , num2str(neuronCoords.x(i)) ]);
                %                  display(['__y:' , num2str(neuronCoords.y(i)) ]);
                %                  display(['__z:' , num2str(neuronCoords.z(i)) ]);
                cylPresent = find(cellfun('length',regexp(neuronPartName,'cyl')) == 1);
                if (cylPresent)
                    neuronCoords.r(i)=eval(propertiesNerve.r);
                    neuronCoords.h(i)=eval(propertiesNerve.h);
                    %                      display(['__r:' , num2str(neuronCoords.r(i)) ]);
                    %                      display(['__h:' , num2str(neuronCoords.h(i)) ]);
                end %if
                
                sphPresent = find(cellfun('length',regexp(neuronPartName,'sph')) == 1);
                if (sphPresent)
                    neuronCoords.r(i)=eval(propertiesNerve.r);
                    neuronCoords.h(i)=eval(propertiesNerve.r);
                    %                      display(['__r:' , num2str(neuronCoords.r(i)) ]);
                end %if
                
                
            end
            
            xNeuronStart = max(neuronCoords.x);  yNeuronStart = max( neuronCoords.y); zNeuronStart = max(neuronCoords.z);
            radiusNerve = neuronCoords.r; heightNerve=sum( neuronCoords.h);
            geometryProperties.nerveLength=heightNerve;
            
            % extract the electrode geometry
            numberOfElectrods=size(electrodeGeometryNames);
            numberOfElectrods=numberOfElectrods(1,1);
            electrodeCoords.x=zeros(numberOfElectrods,1);
            electrodeCoords.y=zeros(numberOfElectrods,1);
            electrodeCoords.z=zeros(numberOfElectrods,1);
            electrodeCoords.radius=zeros(numberOfElectrods,1);
            electrodeCoords.height=zeros(numberOfElectrods,1);
            %electrodeCoords.axis=zeros(numberOfElectrods,1);
            
            for i=1:numberOfElectrods
                electrodeName=electrodeGeometryNames(i,:);
                propertiesElectrode=mphgetproperties(retinaModel.geom('geom1').feature(electrodeName));
                electrodeCoords.x(i)=eval(propertiesElectrode.x);
                electrodeCoords.y(i)=eval(propertiesElectrode.y);
                electrodeCoords.z(i)=eval(propertiesElectrode.z);
                %electrodeCoords.axis(i)=propertiesElectrode.axis;
                
                if (cell2mat(strfind(electrodeName,'cyl')))
                    electrodeCoords.radius(i)=eval(propertiesElectrode.r);
                    electrodeCoords.height(i)=eval(propertiesElectrode.h);
                end %if
            end
            
            geometryProperties.xNeuronStart=xNeuronStart;
            geometryProperties.yNeuronStart=yNeuronStart;
            geometryProperties.zNeuronStart=zNeuronStart;
            geometryProperties.radiusNerve=radiusNerve;
            geometryProperties.electrodeCoords=electrodeCoords;
            geometryProperties.heightNerve=heightNerve;
            geometryProperties.neuronCoords=neuronCoords;
            
            
            display(['init_x coord: ',num2str(geometryProperties.xNeuronStart)]);
            display(['init_y coord: ',num2str(geometryProperties.yNeuronStart)]);
            
            
            
            
        end
        
        
        function geometryProperties=getExVoltages(geometryProperties)
            xNeuronStart = geometryProperties.xNeuronStart;
            yNeuronStart = geometryProperties.yNeuronStart;
            zNeuronStart = geometryProperties.zNeuronStart;
            radiusNerve =  geometryProperties.radiusNerve;
            heightNerve = geometryProperties.heightNerve;
            electrodeCoords = geometryProperties.electrodeCoords;
            neuronCoords = geometryProperties.neuronCoords;
            retinaModel=geometryProperties.model;
            
            
            %% Step 4:
            % find out how the geometry is arranged (how many parallel nerves) %
            nrOfNerveStrains=1;
            compartmentNumbers=100;
            
            %% Step 5:
            % Calculate the pointes of measurement
            geometryProperties.electrodeDistance=abs(zNeuronStart-min(electrodeCoords.z));
            compartmentPerStrain=compartmentNumbers/nrOfNerveStrains;
            delta=heightNerve/(compartmentNumbers-1);
            measurementPoints=zeros(compartmentPerStrain*nrOfNerveStrains,4);
            
            thisZ=zNeuronStart;
            for i=1:length(measurementPoints)-1
                thisX = xNeuronStart;
                thisY = yNeuronStart+(i-1)*delta;
                thisZ = findZinGeometry(thisX,thisY, neuronCoords,thisZ);
                
                measurementPoints(i,:)=[1,thisX,thisY,thisZ];
                %   measurementPoints(i,:)=[1,xNeuronStart,min(neuronCoords.y)+(i-1)*delta,zNeuronStart];
            end %for
            measurementPoints(compartmentNumbers,:)=[1,xNeuronStart,min(neuronCoords.y)+(compartmentNumbers)*delta,zNeuronStart];
            
            % TODO: use the information in the first row, do not simply cut it off!%
            measurementPoints(:,1)=[];
            
            
            %% Step 6
            %get the corresponding values from the FEM model
            
            %Syntax: 3xn => 1xn
            % x1 x2 ...
            % y1 y2 ... => V1 v2
            % z1 z2 ...
            voltageArray=mphinterp(retinaModel,'V','coord',measurementPoints');
            if  (isempty (voltageArray))
                error='can not load external potential. Abort.';
                disp(error);
                return;
            end
            
            geometryProperties.distanceBetrweenCompartments=delta;
            voltageArray=voltageArray*geometryProperties.potMultiplicator;
            geometryProperties.externalVoltages=voltageArray *1000; %get mV
            
            %% Step 7: identify feature for each compartment
            geometryProperties.compFeatures=createCompNeuronTypeVector(geometryProperties,measurementPoints,neuronCoords);
            geometryProperties.compartmentRadius=createNeuronRadiusVector(geometryProperties,radiusNerve);
            
            %% dirty workaround, b/c radius identification does not really work
            geometryProperties.compFeatures = [ repmat({'Soma'},1,1) ; repmat({'IS'},3,1); repmat({'SOCB'},3,1); repmat({'Axon'},93,1) ];
            geometryProperties.compartmentRadius=[ repmat(10*10^-6,1,1) ; repmat(10^-6,3,1); repmat(10^-6,3,1); repmat(10^-6,93,1) ];
            
            
        end % function getExVoltages
        
        
        %% Chose and start a calculation Model (e.g. hodgkin huxley)
        function thismodel=envokeSimulationModel(simulationModel, geometryProperties)
            if (~exist(simulationModel, 'file'))
                error='Model not implemented';
                disp(error);
                return;
            end
            
            % Create a model-Handler of parametered simulation-model
            modelFunctionhandler = str2func(simulationModel);
            thismodel=modelFunctionhandler(geometryProperties);
            
            % Check if model is an isa (InStAnce of) the abstract model
            if (~isa(thismodel,'basicModel'))
                error='No model loaded. exiting';
                disp(error);
                return;
            else
                display(['Using ',simulationModel, ' model']);
            end %if (~isa(thismodel,'abstractModel'))
            
        end % function envokeSimulationModel
        
        %% display all available physical imput geometries on console
        function displayPhysicModels()
            
            listOfImplementedModelsStruct=dir('src/retinaElectrodeModels/geometry*.m');
            for i=1:length(listOfImplementedModelsStruct)
                listOfImplementedModelsStruct(i).geoName=strrep(listOfImplementedModelsStruct(i).name,'geometry_','');
                listOfImplementedModelsStruct(i).geoName=strrep(listOfImplementedModelsStruct(i).geoName,'.m','');
                display(listOfImplementedModelsStruct(i).geoName);
            end
            display(' ');
            display('or in addition, the following test parameters are available: ');
            display('testValues');
            display('testPointStimulation');
            display('NoStimulation');
            display('somaStimulationOnly');
            
            display(' ');
        end %function displayPhysicModels()
        
        
        %% find coords under first electrode
        function [xMin,yMin]=findStartCoords(geometryProperties)
            xMin = min(geometryProperties.electrodeCoords.x)-geometryProperties.electrodeCoords.radius(1)-geometryProperties.mapResolution/1000/1000;
            yMin = min(geometryProperties.electrodeCoords.y)-geometryProperties.electrodeCoords.radius(1)*2;
        end
        
        %% find coords under last electrode
        function [xMax,yMax]=findEndCoords(geometryProperties)
            xMax = geometryProperties.electrodeCoords.x(1)+geometryProperties.mapResolution/1000/1000/2; % only calculate till symmetry axis
            yMax = max(geometryProperties.electrodeCoords.y)+geometryProperties.electrodeCoords.radius(1)*2; % double electrode radius, to make enough space up in the image (expensive!!)
        end
        
        function geometryProperties=adoptNeuronCoords(geometryProperties,xMin,yMin,xMax,yMax,iteration)
            thisResolution=geometryProperties.mapResolution*10^-6; % um => m
            tmp_xNeuronStart = geometryProperties.xNeuronStart;
            tmp_yNeuronStart = geometryProperties.yNeuronStart;
            
            if(iteration==0)
                geometryProperties.xNeuronStart=xMin;
                geometryProperties.yNeuronStart=yMin;
            else
                if(tmp_xNeuronStart+thisResolution>xMax)
                    geometryProperties.xNeuronStart=xMin;
                    
                    if(tmp_yNeuronStart+thisResolution>yMax)
                        geometryProperties.yNeuronStart=yMin;
                    else
                        geometryProperties.yNeuronStart = geometryProperties.yNeuronStart + thisResolution;
                    end %  if(tmp_yNeuronStart+thisResolution>yMax)
                    
                    
                else
                    geometryProperties.xNeuronStart = geometryProperties.xNeuronStart + thisResolution;
                end % if(tmp_xNeuronStart+thisResolution>xMax)
                
                
                
            end % iteration==0)
            display(['new messurement points:']);
            display(['x: ',num2str(geometryProperties.xNeuronStart)]);
            display(['y: ',num2str(geometryProperties.yNeuronStart)]);
            
        end % function adoptNeuronCoords
        
        
        function APfound=findAPinData(Y)
            APfound=false;
            maxPot=max(max(Y(:,1,:)));
            minPot=min(min(Y(:,1,:)));
            baseLine=Y(1,1,1);
            minThreshold=40;
            APfound=~(maxPot<(baseLine+minThreshold));
        end
        
        function done=calculationComplete(geometryProperties,xMin,yMin,iteration)
            done=false;
            tmp_xNeuronStart = geometryProperties.xNeuronStart;
            tmp_yNeuronStart = geometryProperties.yNeuronStart;
            
            if(iteration>0)
                if(tmp_xNeuronStart== xMin && tmp_yNeuronStart==yMin)
                    done=true;
                end
            end % if(iteration!=1)
            
        end
        
        function plotCheck(Y,T,x,y,v)
            return;
            f = figure();
            axisX.min=min(T);
            axisY.min=min(min(Y(:,1,1)))*1.2;
            axisX.max=max(T);
            axisY.max=max(max(Y(:,1,:)))*1.1;
            if axisY.max < 0
                axisY.max=0;
            end
            
            hold on
            compartment=1;
            plot(T,Y(:,1,compartment),'b-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'r-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'g-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'c-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'b.-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'r.-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'g.-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'c.-'); compartment=compartment+10;
            plot(T,Y(:,1,compartment),'y-');
            axis([axisX.min axisX.max axisY.min axisY.max]);
            hold off
            mkdir ('/tmp/mapPlot')
            saveas(f, ['/tmp/mapPlot/mapPlot_x',num2str(x),'_y',num2str(y),'_x',num2str(v),'.png']);
            close all;
        end
        
        
    end %methods (Static)
end %classdef
