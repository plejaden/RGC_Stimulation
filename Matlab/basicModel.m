%% 2016-05-07 13:27:58 
classdef basicModel
    %% Abstract model for definition of MembraneModels
    
    properties (Abstract)
        iPause;
        iStim; %(uA/cm^2) amount of current stimulation for iStimTime
        iStimTime; %(ms) time of current stimulation
        simulationTime; %(s)
        startParameters;
        numberOfStateVariables;
        Ve; r; area; volume;
        ohmAx;
        
        %performance analysis
        t0;
        timeout;
    end % properties (Abstract)
    
    properties (Abstract = false)
        % define the properties of the class here, (like fields of a struct)
        
    end %(Abstract = false)
    
    methods (Abstract=true)
        
        
        
        [result] = h0(this);
        [result] = m0(this);
        [result] = n0(this);
        
        [result] = differentialEq(this,t,y);
        
        [T, Y] = postProcessing(this,T,Y);
        
    end %methods (Abstract=true)
    
    %%
    methods (Abstract = false)
        function [T, Y] = run(this,geometryProperties)
            clear 'currentObserver'; global currentObserver;
            
            % timing analysis
            this.t0=clock;
            
            %% current analysis
            currentObserver.compartment=1; currentObserver.time = 0; currentObserver.INa = 0; currentObserver.IK = 0; currentObserver.externalVoltage=geometryProperties.externalVoltages;
            currentObserver.IKA = 0; currentObserver.IKCa = 0; currentObserver.ICa = 0; currentObserver.IL = 0; currentObserver.stimulusCurrent = 0; currentObserver.externalCurrentSet=false;
            
            % rotate the external voltages, if they are not in array shape
            voltageForm=size(geometryProperties.externalVoltages);
            if (voltageForm(1)<voltageForm(2))
                geometryProperties.externalVoltages=geometryProperties.externalVoltages';
            end
            global activation;
            activation.current=zeros(max(voltageForm),1);
            
            % rho = R * q / l => R = L * rho / q
            rhoFiber = this.rhoFiber*10^-3; %[Ohm*mm]
            compLength=geometryProperties.distanceBetrweenCompartments; % [m]
            compRad=geometryProperties.compartmentRadius; % [m]
            compNumb=geometryProperties.numberOfCompartments;
            thisSurface=2.*compRad.*pi().*compLength+2.*compRad.*compRad.*pi();
            display(['compLength:', num2str(compLength)]);
            display(['compSurface (first cmp): ', num2str(thisSurface(1)), ' m^2']);
            display(['compSurface (first cmp): ', num2str(thisSurface(1)*(1000*1000)), ' mm^2']);
            display(['compSurface (first cmp): ', num2str(thisSurface(1)*(1000*1000*1000*1000)), ' um^2']);
            
            l_a=(compLength * 2); % [m]
            A_a=(compRad .* 1000).^2*pi; % [mm^2]
            
            this.ohmAx= rhoFiber*l_a./A_a; % ( [Ohm*mm] ) * m / (mm^2) = [kOhm]
            
            
            % save the geometry Properties to external file: very useful for debugging %
            geometryProperties.modelproperties=this;
            
            simTime=[0 this.simulationTime];
            
            try
                %% solve the equation
                % preferred solver: ode113
                %            [T, Y] = ode15s(@(t,y) differentialEq(this, t, y ), simTime, this.startParameters); % stiff
                %            [T, Y] = ode23(@(t,y) differentialEq(this, t, y ), simTime, this.startParameters); % non stiff
                %            [T, Y] = ode23s(@(t,y) differentialEq(this, t, y ), simTime, this.startParameters); % stiff
                %            [T, Y] = ode45(@(t,y) differentialEq(this, t, y ), simTime, this.startParameters);% non stiff
                [T, Y] = ode113(@(t,y) differentialEq(this, t, y ), simTime, this.startParameters);% stiff
            catch
                display('Error: No solution detected');
                global exceptionValue;
                
                Y=exceptionValue;
                T=zeros(compNumb);
                
            end
            
            %%
            if (~isreal(Y))
                disp('Y is complex, this should not happen!!!');
                return;
            end
            
            if (size(Y,2) <=compNumb)
                disp('Y does not contain gating information, this should not happen!!!');
                return;
            end
            
            % Y is an array (Timesteps x (compartments * variables) e.g. 4612 x 88 %
            % But we want: Timesteps x variables x compartments e.g. 4612 x 8 x 11 %
            calcPointsPerTime=max(size(Y(1,:))); % hols: |comparment| * |statvariables|%
            compL=calcPointsPerTime/(this.numberOfStateVariables+1);
            timeSlices=length(T);
            tempY=zeros(timeSlices,this.numberOfStateVariables+1, compL);
            
            for i=1:this.numberOfStateVariables+1;
                tempY(:,i,:)=Y(:,(i-1)*compL+1:i*compL);
            end % for i=1:this.numberOfStateVariables+1;
            Y=tempY;
            
            % do some post-processing
            [T, Y]=this.postProcessing(this,T,Y);
            
        end %function [T, Y] = run(this)
    end %methods (Abstract = false)
    
    methods (Static)
        
        %% numberOfCompartments [1]; Ve [mV]; pulseMagnification [1]
        function stimulusCurrentExtracomp = compartmentStim(this,numberOfCompartments, Ve, pulseMagnification)
            global currentObserver
            stimulusCurrentExtracomp=zeros(numberOfCompartments,1);
            if (pulseMagnification~=0)
                
                %first compartment
                if (this.iStim==0 )
                    deltaVEfirst=Ve(1)-Ve(2); % mV
                    firstOhmAx=(this.ohmAx(1)+this.ohmAx(2))./2;
                    stimulusCurrentExtracomp(1) = deltaVEfirst./(firstOhmAx.*2);  % mV / kOhm = uA
                end % (iStim~=0)
                
                % middle compartments
                deltaVE1=Ve(1:numberOfCompartments-2)-Ve(2:numberOfCompartments-1);
                deltaVE2=Ve(3:numberOfCompartments)-Ve(2:numberOfCompartments-1); % the delta values, shifted at once.
                centerOhmAx=(this.ohmAx(1:numberOfCompartments-2)+this.ohmAx(2:numberOfCompartments-1))/2';
                stimulusCurrentExtracomp(2:numberOfCompartments-1) = (deltaVE1 + deltaVE2)./(centerOhmAx');  % mV / kOhm = uA ; the vector is rotated automatically %
                
                % last compartment
                deltaVElast=Ve(numberOfCompartments-1)-Ve(numberOfCompartments);
                lastOhmAx=this.ohmAx(numberOfCompartments-1)+this.ohmAx(numberOfCompartments);
                stimulusCurrentExtracomp(numberOfCompartments) = deltaVElast./(lastOhmAx);  % mV / kOhm = uA
                stimulusCurrentExtracomp = stimulusCurrentExtracomp.*pulseMagnification;
                
                if (currentObserver.externalCurrentSet==false)
                    currentObserver.externalVoltage=Ve;
                    currentObserver.externalCurrent=stimulusCurrentExtracomp;
                    currentObserver.externalCurrentSet=true;
                end
                
            end % if (pulseMagnification~=0)
            
        end % function compartmentStim
        
        
        function X = rStim(pause,current)
            if (pause~=0)
                X=current;
            else
                X=0;
            end
            
        end % function rStim
        
        %% current calculations
        function stimulusCurrent = linkCompartments(this,compL,time,V,stimulationDesc,iStim,Ve)
            
            global activation;
            stimulusCurrentIntercomp=zeros(compL,1);
            stimulusCurrentExtracomp=zeros(compL,1);
            
            pulsTime1 = stimulationDesc(1,1);           pulseMagnif1 = stimulationDesc(1,2);
            pulsTime2 = stimulationDesc(2,1)+pulsTime1; pulseMagnif2 = stimulationDesc(2,2);
            pulsTime3 = stimulationDesc(3,1)+pulsTime2; pulseMagnif3 = stimulationDesc(3,2);
            pulsTime4 = stimulationDesc(4,1)+pulsTime3; pulseMagnif4 = stimulationDesc(4,2);
            
            %% Apply somatic current, if enabled
            % The first compartment gets the stimulus from Istim, all other
            % compartments get the stimulus from their neighbor comp %
            %
            % istim->[Soma][][][][][]
            %
            if (iStim~=0)
                % stimulusCurrentIntercomp(1)=this.rStim(time,pulsTime1, pulsTime2, iStim);
                if (time < pulsTime1)
                    %stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif1);
                    stimulusCurrentIntercomp(1)=pulseMagnif1.*iStim;
                elseif (time < pulsTime2)
                    %stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif2);
                    stimulusCurrentIntercomp(1)=pulseMagnif2.*iStim;
                elseif (time < pulsTime3)
                    %stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif3);
                    stimulusCurrentIntercomp(1)=pulseMagnif3.*iStim;
                elseif (time < pulsTime4)
                    %stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif4);
                    stimulusCurrentIntercomp(1)=pulseMagnif4.*iStim;
                elseif (time > pulsTime4)
                    %stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, 0);
                    stimulusCurrentIntercomp(1)=0;
                end % if (time < pulsTime1)
                
            else
                stimulusCurrentIntercomp(1)=0;
                deltaVfirst=-V(1)+V(2);
                lastOhmAx=(this.ohmAx(compL-1)+this.ohmAx(compL))/2;
                stimulusCurrentIntercomp(1) = deltaVfirst./(lastOhmAx); % mV / kOhm = ?A
                stimulusCurrentIntercomp(1) = stimulusCurrentIntercomp(1)/2;
            end % (iStim~=0)
            
            %% link the compartments
            % calculate the current between compartments
            % [][mV]-[mV][][]
            
            if (compL>1)
                % middle compartments
                deltaV1=V(1:compL-2)-V(2:compL-1);
                deltaV2=V(3:compL)-V(2:compL-1); % the delta values, shifted at once.
                centerOhmAx=(this.ohmAx(1:compL-2)+this.ohmAx(2:compL-1))/2;
                stimulusCurrentIntercomp(2:compL-1) = (deltaV1 + deltaV2)./(centerOhmAx.*2); % mV / kOhm
                
                %last compartment
                deltaVlast=V(compL-1)-V(compL);
                lastOhmAx=(this.ohmAx(compL-1)+this.ohmAx(compL))/2;
                stimulusCurrentIntercomp(compL) = deltaVlast/(lastOhmAx); % mV / kOhm
                
            end % if (compL>1)
            
            %% Apply the external activating current
            % Use the Ve for getting the influence of the external voltage on the core of the compartment
            %
            % ----Ve---
            % [ ][Vi][]
            if (compL>1)
                if (time < pulsTime1)
                    stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif1);
                elseif (time < pulsTime2)
                    stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif2);
                elseif (time < pulsTime3)
                    stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif3);
                elseif (time < pulsTime4)
                    stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, pulseMagnif4);
                elseif (time > pulsTime4)
                    stimulusCurrentExtracomp = this.compartmentStim(this,compL,Ve, 0);
                end
            end % if (compL>1)
            %stimulusCurrentExtracomp=0;
            
            %% Sum up the results
            stimulusCurrent=stimulusCurrentExtracomp+stimulusCurrentIntercomp;
            
            if (activation.current==0)
                activation.current=stimulusCurrent;
                if (activation.current~=0) display(['Current measured at ms: ' num2str(time) ]); end
            end % if (activation.current==0)
        end % function linkCompartments
        
        function X = surfaceToVolumeRatio(r,l)
            thisVolume=r.*r*pi()*l;
            surface=2.*r.*pi().*l+2.*r.*r.*pi();
            X=surface./thisVolume;
        end %function
        
    end % methods (Static)
end % classdef basicModel
