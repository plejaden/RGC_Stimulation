classdef noModel < basicModel
    properties
        simulationTime=10;
        temperatureCelsius=20;
        
        t0; timeout=300; %s
        
        C = 1.0; %(uF/cm^2) Membrane capacitance
        Vrest=-70; %[mV]
        VNa = 115; %(mV) Sodium Voltage
        VK = -12; %(mV) Potassium Voltage
        VL = 10.613; %(mV) Leakage Voltage
        gNa = 120; %(mS/cm^2) Sodium conductance
        gL = 0.3; %(mS/cm^2) Leakage conductance
        gK = 36; %(mS/cm^2) Potassium conductance
        
        %% simulation
        iPause=0.5;
        iStim; %(uA/cm^2) Amount of current stimulation for 0.12msec
        iStimTime=0.120;
        
        % biphasic stimulation: for this time, the mentioned multiplyer will be applied %
        % format: (time in [ms], pot multiplicator [1] ; )
        % [ puls1,pot1 ; puls2,pot2  ; puls3,pot3 ; puls4,pot4]
        % exmaple: [0.5,0.15,0.2,0.1; 0,-1,0,+1]
        
        % intracellular stimulation:
        stimulationDescInt=[0.5,0; 0.2,1; 3.5,0; 0.2,0];
        
        % extracellular stimulation:
        stimulationDescExt=[0.5,0; 0.2,1; 1,0; 0.4,-0.5]; % with electrode uncharge
        % stimulationDescExt=[0.5,0; 0.2,1; 1,0; 0.4,0]; % without electrode uncharge
        
        startParameters;
        numberOfStateVariables=3; % h m n
        Ve=0;  r; area; volume; k_q10;
        rhoFiber = 1.1; %[Ohm*m]
        stimulationDesc; ohmAx;
        
        T0=6.3; % ?C
        q10=3;
        
        
        
    end %properties
    
    methods
        
        % The standard constructor sets the temperature of the model.
        function this = noModel(geometryProperties)
           
            this.Ve=geometryProperties.externalVoltages;
            this.iStim=geometryProperties.iStim;
            geometryProperties.compartmentRadius=geometryProperties.compartmentRadius/10; % radius in cm
            this.r=geometryProperties.compartmentRadius;
            this.area=geometryProperties.compartmentRadius.^2*pi;
            this.volume=geometryProperties.distanceBetrweenCompartments*this.area;
            
            if (this.iStim==0)
                display('using extracellular stimulation function');
                this.stimulationDesc=this.stimulationDescExt;
            else
                display('using intercellular stimulation function');
                this.stimulationDesc=this.stimulationDescInt;
            end % (this.iStim~=0)
            
            % temperature factor
            this.k_q10 = this.q10^((this.temperatureCelsius-this.T0)/10);
            display(['k_q10:', num2str(this.k_q10)]);
            
            l=geometryProperties.numberOfCompartments; % get the number of compartments from the voltage array%
                        
            V0=repmat(this.Vrest,1,l)';
            this.startParameters=[V0 this.h0(l) this.m0(l) this.n0(l)];
            %obj.k=temp-263.15;
        end
        
        function X = differentialEq(this,time,x)
            global currentObserver;
            compL=length(x)/(this.numberOfStateVariables+1);
            
            %% allocate memory
            X=zeros(compL,4);
            
            currentObserver.time = [currentObserver.time; time]; currentObserver.INa = [currentObserver.INa; 0];
            currentObserver.IK = [currentObserver.IK; 0]; currentObserver.IL = [currentObserver.IL; 0];
            
            X=X(:);
        end
        
        function X = h0(this,l)
            X=this.alpha_h(0)/(this.alpha_h(0)+this.beta_h(0));
            X=repmat(X,1,l)';
        end %function
        
        function X = m0(this,l)
            X=this.alpha_m(0)/(this.alpha_m(0)+this.beta_m(0));
            X=repmat(X,1,l)';
        end %function
        
        function X = n0(this,l)
            X=this.alpha_n(0)/(this.alpha_n(0)+this.beta_n(0));
            X=repmat(X,1,l)';
        end %function
        
        function X = p0(~,~)
            X=0; %% not needed in HH
        end %function
        
    end %methods
    
    methods (Static)
        
        function X = E(this,V) %% todo: check if necessary
            X=V - this.Vrest;
        end %function E
        
        function X = alpha_h(V)
            X=0.07.*exp(-V./20);
        end %function
        
        function X = beta_h(V)
            X=1./(exp(0.1.*(29.99-V))+1);
        end %function
        
        function X = alpha_m(V)
            X=(0.1.*(25-V))./(exp(0.1.*(25-V))-1);
        end %function
        
        function X = beta_m(V)
            X=4.*exp(-V/18);
        end %function
        
        function X = alpha_n(V)
            X=(0.1.*0.1.*(10-V))./(exp(0.1.*(10-V))-1);
        end %function
        
        function X = beta_n(V)
            X=0.125.*exp(-V/80);
        end %function
        
        function X = alpha_p(V)
            X=V*0; %% Not needed in HH
        end %function
        
        function X = beta_p(V)
            X=V*0; %% Not needed in HH
        end %function
        
        function [T, Y] = postProcessing(this,T,Y)
            return;
            plot(T,Y(:,1,1));
            xlabel('t [ms]'); ylabel('membran potential [mV]');
            filename=[ 'hh_12degree' '+' num2str(this.iStim) 'ua+' num2str(this.stimulationDesc(3,1)) 'ms+' datestr(now, 'yyyy-mm-dd_HH:MM:SS')  ];
            mkdir('/tmp/membranPotential');
            myStyle = hgexport('factorystyle'); myStyle.Format = 'png'; myStyle.Width = 3; myStyle.Height = 1; myStyle.Resolution = 300; myStyle.Units = 'inch'; myStyle.FixedFontSize = 12;
            hgexport(gcf, [ '/tmp' '/' 'membranPotential' '/' filename '.png' ], myStyle, 'Format', 'png');
            close all;
        end %function postProcessing
        
    end
end %classdef
