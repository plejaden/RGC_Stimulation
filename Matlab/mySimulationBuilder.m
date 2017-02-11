classdef mySimulationBuilder
   % Author: Watzinger Anton
   % Last edited: 2017-02-11 
    
    %% Start with:
    % clear; mySimulationBuilder('fcn2010_cat_20deg',1,'2electrodes_epiretinal',100);
    methods
        % class constructor
        function this = mySimulationBuilder(simulationModelName, potentialMultiplicator, electrodeConfiguration,compartmentNumbers)
            
            cd /home/plejaden/Dokumente/Uni/MasterThesis/trunk;
            
            global error;
            error = ''; % no error occured
            starttime=datestr(now, 'yyyy-mm-dd HH:MM:SS');
            display(' ');
            display([ 'Simulation started at ' starttime]);
            close all;
            warning ('off','all');
            
            tic;
            numberOfParameters=nargin();
            
            if(numberOfParameters <= 3)
                display('Quick-guide for using this system');
                display('mySimulationBuilder(''fcn2010_cat'',1,''2electrodes_epiretinal'',100)');
                display(' ');
                display('#1 parameter: used model for calculation. Implemented models are:');
                dir('membraneModels/*.m')
                display(' ');
                display('#2 parameter: voltage multiplier. use 1 to use the exact volatge from the comsol model');
                display(' ');
                display('#3 parameter: Geometry. Implemented geometries are:');
                this.displayPhysicModels();
                display(' ');
                display('#4 parameter: Number of compartments. use a value from 10 to 70');
                return;
            end % 
            
            
            display('available Geometries:');
            this.displayPhysicModels();
            
            if(numberOfParameters == 4)
                
                % create retina, nerve and electrode geometry, mesh, run physics %
                if (isempty (error))
                    display('## calculate finite element potentials ##'); toc;
                    geometryProperties=this.organizeModel(this,electrodeConfiguration,compartmentNumbers,potentialMultiplicator);
                    geometryProperties.filename=[ simulationModelName,'+', num2str(potentialMultiplicator),'+', electrodeConfiguration,'+',num2str(compartmentNumbers),'+',datestr(now, 'yyyy-mm-dd_HH:MM:SS') ];
                    %display(geometryProperties.externalVoltages);
                end
                
                % envokes the nerve model for calculating the membrane
                % potential
                if (isempty (error))
                    simulationModel = this.envokeSimulationModel(simulationModelName, geometryProperties);
                end %  if (isempty (error))
                
                %Calculate the model
                if (isempty (error))
                    display(' '); display('## solveing ode ##'); toc;
                    [T, Y] = simulationModel.run(geometryProperties);
                end % if (isempty (error))
                
                if (isempty (error))
                    % check if there are some differences in the resulting matrix, which may be an AP %
                    if (max(max(Y(:,1,:)))<min(min(Y(:,1,:)))+50)
                        display('WARNING: No AP spike detected!'); toc;
                    end % if (max(max(Y(:,1,:)))<min(min(Y(:,1,:)))+50)
                end % if (isempty (error))
                
                %Plot the results
                global currentObserver;
                filename=['output/',geometryProperties.filename,'.mat'];
                save (filename);
                sliderPlot(T,Y,geometryProperties);
                
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
                
            end % if strcmpi(modelType,'noStimulation')
            
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

            end % if strcmpi(modelType,'intercellularStimulationOnly')
            
            
            % no external stimulation
            if strcmpi(modelType,'intercellularSinglecompStimulationOnly')
                display('_-_-_-_-_-_-')
                display('INTERCELLULAR-STIMULATION (SINGLE COMP):');
                display('Disabled all stimulus but the "potentialMultiplicator" [ua/cm^2] at the first compartment');
                geometryProperties.model='noModel-virtual soma stimulaton only';
                geometryProperties.rho = 50.5;
                geometryProperties.numberOfCompartments=1;
                geometryProperties.iStim=potentialMultiplicator;
                
                geometryProperties.compFeatures =  {'Soma'} ;
                geometryProperties.compartmentRadius= (5*10^-6) ;
                geometryProperties.electrodeDistance=6*10^-5; % not used
                geometryProperties.distanceBetrweenCompartments=30*10^-6; %m
                
                geometryProperties.externalVoltages=zeros(1,1);
                return;
            end %   if strcmpi(modelType,'intercellularSinglecompStimulationOnly')
            
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
            end % if (exist(geoFuncName, 'file') && exist(gdescName, 'file') )
            
            % envoking the gmodel and gdesc functionpointers.
            geometryProperties.model=gmodel();
            geometryProperties.desc=gdesc(); toc;
            
            geometryProperties.numberOfCompartments=compartmentNumbers;
            geometryProperties=this.getExVoltages(geometryProperties); % add most of the properties
        end % geometryProperties=organizeModel(this,modelType,compartmentNumbers,potentialMultiplicator)
        
        %% get the external Voltages form the geometry
        function geometryProperties=getExVoltages(geometryProperties)
            
            retinaModel=geometryProperties.model;
            compartmentNumbers=geometryProperties.numberOfCompartments;
            geometryProperties.distanceBetrweenCompartments=4; % ??
            
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
                display(['reading neuronal part: ',neuronPartName]);
                neuronCoords.name{i}=neuronPartName;
                neuronCoords.x(i)=eval(propertiesNerve.x);
                neuronCoords.y(i)=eval(propertiesNerve.y);
                neuronCoords.z(i)=eval(propertiesNerve.z);
                
                cylPresent = find(cellfun('length',regexp(neuronPartName,'cyl')) == 1);
                if (cylPresent)
                    neuronCoords.r(i)=eval(propertiesNerve.r);
                    neuronCoords.h(i)=eval(propertiesNerve.h);
                end %if (cylPresent)
                
                sphPresent = find(cellfun('length',regexp(neuronPartName,'sph')) == 1);
                if (sphPresent)
                    neuronCoords.r(i)=eval(propertiesNerve.r);
                    neuronCoords.h(i)=eval(propertiesNerve.r);
                end %if (sphPresent)
                
            end % for i=1:numberOfNeuronparts
            
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
            
            for i=1:numberOfElectrods
                electrodeName=electrodeGeometryNames(i,:);
                propertiesElectrode=mphgetproperties(retinaModel.geom('geom1').feature(electrodeName));
                electrodeCoords.x(i)=eval(propertiesElectrode.x);
                electrodeCoords.y(i)=eval(propertiesElectrode.y);
                electrodeCoords.z(i)=eval(propertiesElectrode.z);
                if (cell2mat(strfind(electrodeName,'cyl')))
                    electrodeCoords.radius(i)=eval(propertiesElectrode.r);
                    electrodeCoords.height(i)=eval(propertiesElectrode.h);
                end %(cell2mat(strfind(electrodeName,'cyl')))
            end % i=1:numberOfElectrods
            
            %% Step 4:
            % find out how the geometry is arranged (how many parallel nerves) %
            nrOfNerveStrains=1;
            
            %% Step 5:
            % Calculate the pointes of measurement
            geometryProperties.electrodeDistance=abs(zNeuronStart-min(electrodeCoords.z));
            compartmentPerStrain=compartmentNumbers/nrOfNerveStrains;
            delta=heightNerve/(compartmentNumbers-1);
            measurementPoints=zeros(compartmentPerStrain*nrOfNerveStrains,4);
            
            thisZ=zNeuronStart;
            for i=1:length(measurementPoints)-1
                thisX = xNeuronStart;
                thisY = min(neuronCoords.y)+(i-1)*delta;
                thisZ = findZinGeometry(thisX,thisY, neuronCoords,thisZ);
                
                measurementPoints(i,:)=[1,thisX,thisY,thisZ];
            end %for i=1:length(measurementPoints)-1
            measurementPoints(compartmentNumbers,:)=[1,xNeuronStart,min(neuronCoords.y)+(compartmentNumbers)*delta,zNeuronStart];
            
            % multi-neuron in row 1, but currently not supported.
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
            end % (isempty (voltageArray))
            
            geometryProperties.distanceBetrweenCompartments=delta;
            voltageArray=voltageArray*geometryProperties.potMultiplicator;
            geometryProperties.externalVoltages=voltageArray *1000; %V => mV
            
            %% Step 7: identify feature for each compartment
            geometryProperties.compFeatures=createCompNeuronTypeVector(geometryProperties,measurementPoints,neuronCoords);
            geometryProperties.compartmentRadius=createNeuronRadiusVector(geometryProperties,radiusNerve);
            
        end
        
        %% Chose and start a membrane model (e.g. hodgkin huxley)
        function thismodel=envokeSimulationModel(simulationModel, geometryProperties)
            if (~exist(simulationModel, 'file'))
                error='Model not implemented';
                disp(error);
                return;
            end % if (~exist(simulationModel, 'file'))
            
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
        
    end %methods (Static)
end %classdef
