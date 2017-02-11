classdef fcn2010_cat_20deg_singleComp < basicModel
%% 2016-05-07 13:27:33    
    properties
        simulationTime=5; %ms
        t0; timeout=60; %s
        
        gasConstant=8.31441; % (J/mol*K)
        Faraday=96485.3365; %c/mol
        RTdivF; %mV
        area;
        volume; Temp;
        gNa; gK; gCa;
        T0=35; T=23.5;
        
        %% Temperature dependent factors
        factor_Q10_Na_kin; factor_Q10_K_kin; factor_Q10_Ca_kin
        factor_Q10_Na_gBar; factor_Q10_K_gBar; factor_Q10_Ca_gBar;
        
        %% channel densities soma
        gNaSoma=69.4;%(mS/cm^2)
        gKSoma=32;% (mS/cm2)
        gCaSoma=1.39;% (mS/cm2)
        
        %% channel densities initial segment
        gNaIS=100;%(mS/cm^2)
        gKIS=50.10;% (mS/cm2)
        gCaIS=0.836;% (mS/cm2)
        
        %% channel densities axon
        gNaAxon=124;%(mS/cm^2)
        gKAxon=50;% (mS/cm2)
        gCaAxon=04;% (mS/cm2)
        
        gKCa=0.05;% (mS/cm2)
        gL=0.1; % (mS/cm2)
        
        
        %% elibrium potentials
        VNa = 60.60; %(mV)
        VK = -101.34; % (mV)
        VL= -65; %(mV)
        
        %% Calcium properties
        CaRes=10^-7*10^3; % in M/dm3 => mM/cm3
        CaDiss=10^-6*10^3; % in M/dm3 => mM/cm3
        Cao=1.8; % mM
        
        %% misc properties
        Vrest= -70; % (mV)
        Rn = 150*10^6; %(ohm) = 150 MOhm, input resistance
        r;compartmentLength;
        C = 1; %(mF/cm2)
        TauCA=1.5; %
        
        %% simulation
        iPause=0.5; % ms
        iStimTime=0.120;
        % would be; [0.5,0.12,0,0; 0,1,0,0]
        
        % soma stimulation
        iStim=0; %(uA/cm^2) Amount of current stimulation
        
        % biphasic stimulation: for this time, the mentioned multiplyer will be applied %
        % format: (time in [ms], pot multiplicator [1] ; )
        % [ puls1,pot1 ; puls2,pot2  ; puls3,pot3 ; puls4,pot4]
        % exmaple: [0.5,0.15,0.2,0.1; 0,-1,0,+1]
        
        % intracellular stimulation:
        stimulationDescInt=[0.5,0; 0.2,1; 2.25,0; 0.2,1];
        
        % extracellular stimulation:
       % stimulationDescExt=[0.5,0; 0.2,1; 0.5,0; 0.4,-0.5]; % with electrode uncharge
        stimulationDescExt=[0.5,0; 0.2,1; 1,0; 0.4,0]; % without electrode uncharge
        
        startParameters;
        numberOfStateVariables=5; % h m n c Cai
        Ve=0;
        rhoFiber = 1.1; %[Ohm*m]
        stimulationDesc; ohmAx;
     

        
    end %properties
    
    methods
        
        %% The standard constructor sets the temperature of the model.
        function this = fcn2010_cat_20deg_singleComp(geometryProperties)
            
            global currentObserver;
            global temperature;
            
            this.Temp=273.15+this.T;
            temperature = this.T;
            
            %% Temperature dependent factors
            this.factor_Q10_Na_kin=1.95^((this.T-this.T0)/10);
            this.factor_Q10_K_kin=1.9^((this.T-this.T0)/10);
            this.factor_Q10_Ca_kin=1.95^((this.T-this.T0)/10);
            
            this.factor_Q10_Na_gBar=1.64^((this.T-this.T0)/10);
            this.factor_Q10_K_gBar=1.49^((this.T-this.T0)/10);
            this.factor_Q10_Ca_gBar=1.64^((this.T-this.T0)/10);
            
            %% channel densities soma
            this.gNaSoma=this.gNaSoma*this.factor_Q10_Na_gBar;%(mS/cm^2)
            this.gKSoma= this.gKSoma*this.factor_Q10_K_gBar;% (mS/cm2)
            this.gCaSoma=this.gCaSoma*this.factor_Q10_Ca_gBar;% (mS/cm2)
            
            %% channel densities initial segment
            this.gNaIS=this.gNaIS*this.factor_Q10_Na_gBar;%(mS/cm^2)
            this.gKIS=this.gKIS*this.factor_Q10_K_gBar;% (mS/cm2)
            this.gCaIS=this.gCaIS*this.factor_Q10_Ca_gBar;% (mS/cm2)
            
            %% channel densities axon
            this.gNaAxon=this.gNaAxon*this.factor_Q10_Na_gBar;%(mS/cm^2)
            this.gKAxon=this.gKAxon*this.factor_Q10_K_gBar;% (mS/cm2)
            this.gCaAxon=this.gCaAxon*this.factor_Q10_Ca_gBar;% (mS/cm2)
            
            
            this.Ve=geometryProperties.externalVoltages;
            this.iStim=geometryProperties.iStim;
            this.r=geometryProperties.compartmentRadius*100; % m => cm
            this.compartmentLength=geometryProperties.distanceBetrweenCompartments*100; % m => cm
            this.area=this.r.^2*pi;   % area in cm2
            this.volume=(this.compartmentLength)*this.area;
            
            if (this.iStim==0)
                display('using extracellular stimulation function');
                this.stimulationDesc=this.stimulationDescExt;
            else
                display('using intercellular stimulation function');
                this.stimulationDesc=this.stimulationDescInt;
            end % (this.iStim~=0)
            
            % channel densities vector
            this.gNa = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'SOCB' 'Axon'},[this.gNaSoma,this.gNaIS,this.gNaAxon.*10,this.gNaAxon]);
            this.gK = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'SOCB' 'Axon'},[this.gKSoma,this.gKIS,this.gKAxon,this.gKAxon]);
            this.gCa = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'SOCB' 'Axon'},[this.gCaSoma,this.gCaIS,this.gCaAxon,this.gCaAxon]);
            
            % get the number of compartments from the voltage array%
            numberOfCmp=geometryProperties.numberOfCompartments;
            
            %% current analysis
            currentObserver.compartment=1; currentObserver.time = 0; currentObserver.INa = 0; currentObserver.IK = 0;
            currentObserver.IKA = 0; currentObserver.IKCa = 0; currentObserver.ICa = 0; currentObserver.IL = 0;
            
            V0=repmat(this.Vrest,1,numberOfCmp)';
            this.startParameters=[V0 this.h0(numberOfCmp) this.m0(numberOfCmp) this.n0(numberOfCmp) this.c0(numberOfCmp) this.Cai0(numberOfCmp)];
            this.RTdivF=(this.gasConstant*this.Temp)/this.Faraday;
            
        end
        
        %% the differential equations
        function X = differentialEq(this,time,x)
            compL=length(x)/(this.numberOfStateVariables+1);
            
            global currentObserver;
            
            %prevents error if t0 is not initialized
            if (isempty(this.t0)) this.t0=clock; end
            
            ms = round(etime(clock,this.t0) * 1000);
            if(ms>this.timeout*1000)
                display(['Timeout, maxPotential:',num2str(max(x(1:compL)))]);
                global exceptionValue; exceptionValue=x(1:compL);
                return;
            end
            
            % slice up the column vector to a 'compL'x6 array
            
            V=x(1:compL); h=x(compL+1:2*compL); m=x(2*compL+1:3*compL); n=x(3*compL+1:4*compL);
            c=x(4*compL+1:5*compL); Cai=abs(x(5*compL+1:6*compL));
            thisVe=this.Ve;
            
            %% get the stimulus power for each compartment, if the neuron gets stimulated in first compartment
            stimulusCurrent=this.linkCompartments(this,compL,time,V,this.stimulationDesc,this.iStim, thisVe); % current in mA %
            
            %% allocate memory
            X=zeros(compL,6);
            
            %% Main equation
            thisINa=zeros(compL,1); thisIK=zeros(compL,1); thisIKCa=zeros(compL,1);
            thisICa=zeros(compL,1); thisIL=zeros(compL,1);
            
            % to disable single currents, just comment the according line
            thisINa=this.INa(V,m,h);
            thisIK=this.IK(V,n);
            thisIKCa=this.IKCa(V,Cai);
            thisICa=this.ICa(V,c,Cai);
            thisIL=this.IL(V);
            
            currentObserver.time = [currentObserver.time; time]; currentObserver.INa = [currentObserver.INa; thisINa(currentObserver.compartment)];
            currentObserver.IK = [currentObserver.IK; thisIK(currentObserver.compartment)];
            currentObserver.IKCa = [currentObserver.IKCa; thisIKCa(currentObserver.compartment)]; currentObserver.ICa = [currentObserver.ICa; thisICa(currentObserver.compartment)];
            currentObserver.IL = [currentObserver.IL; thisIL(currentObserver.compartment)];
            
            X(:,1)=(-thisINa-thisIK-thisIKCa-thisICa-thisIL+stimulusCurrent)./this.C;%Vm [mV];

            X(:,2)=(-(this.alpha_h(this,V)+this.beta_h(this,V)).*h+this.alpha_h(this,V)).*this.factor_Q10_Na_kin; %h
            X(:,3)=(-(this.alpha_m(this,V)+this.beta_m(this,V)).*m+this.alpha_m(this,V)).*this.factor_Q10_Na_kin; %m
            X(:,4)=(-(this.alpha_n(this,V)+this.beta_n(this,V)).*n+this.alpha_n(this,V)).*this.factor_Q10_K_kin; %n
            X(:,5)=(-(this.alpha_c(this,V)+this.beta_c(this,V)).*c+this.alpha_c(this,V)).*this.factor_Q10_Ca_kin; %c
            X(:,6)=(this.surfaceToVolumeRatio(this.r,this.compartmentLength).*(-this.ICa(V,c,Cai)./(this.Faraday*2)) .* 10^-3)-((Cai-this.CaRes)./this.TauCA);%Cai
            
            %create a column vector
            X=X(:);
        end
        
        %% Helper functions
        function X=INa(this,V,m,h) % check
            X=this.gNa.*m.^3.*h.*(V-this.VNa);
        end %function INa
        
        function X=IK(this,V, n)
            X=this.gK.*n.^4.*(V-this.VK);
        end %function IK
        
        function X=IKCa(this,V,Cai)
            X=this.gKCaBar(this,Cai) .* (V-this.VL);
        end %function IKCa
        
        function X=ICa(this,V,c,Cai) % check
            relative_VCa=this.VCa(Cai);
            X=this.gCa.*c.^3.*(V-relative_VCa);
        end %function ICa
        
        function X=VCa(this, Cai)
            RTdivF2=(this.gasConstant*this.Temp)/(this.Faraday.*2);
            X=RTdivF2.*log(this.Cao./Cai); % in volt
            X=X*10^3; % V => mV
        end %function VCa
        
        % Leak current
        function X=IL(this,V)
            X=this.gL *(V-this.VL);
        end %function IL
        
        %% start Values
        function X = h0(this,l)
            X=this.alpha_h(this,this.Vrest)./(this.alpha_h(this,this.Vrest)+this.beta_h(this,this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = m0(this,l)
            X=this.alpha_m(this,this.Vrest)./(this.alpha_m(this,this.Vrest)+this.beta_m(this,this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = n0(this,l)
            X=this.alpha_n(this,this.Vrest)./(this.alpha_n(this,this.Vrest)+this.beta_n(this,this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = c0(this,l)
            X=this.alpha_c(this,this.Vrest)./(this.alpha_c(this,this.Vrest)+this.beta_c(this,this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = Cai0(~,l)
            X=1.0000e-04;
            X=repmat(X,1,l)';
        end %function
        
    end %methods
    
    methods (Static)
        function X = alpha_h(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=this.factor_Q10_Na_kin*1.869.*exp(-(V+55)/20);
            elseif(temperature>30)
                X=this.factor_Q10_Na_kin*1.817.*exp(-(V+52)/20);
            end
        end %function
        
        function X = beta_h(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=this.factor_Q10_Na_kin*28.04./(1+exp(-0.1*(V+25)));
            elseif(temperature>30)
                X=this.factor_Q10_Na_kin*27.25./(1+exp(-0.1*(V+22)));
            end
        end %function
        
        function X = alpha_m(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=(this.factor_Q10_Na_kin*-2.804*(V+35))./(exp(-0.1*(35+V))-1);
            elseif(temperature>30)
                X=(this.factor_Q10_Na_kin*-2.725*(V+35))./(exp(-0.1*(35+V))-1);
            end
        end %function
        
        function X = beta_m(this,V)
           temperature = this.T;
            X=0;
            if (temperature<24)
                X=this.factor_Q10_Na_kin*93.46*exp(-(V+60)/18);
            elseif(temperature>30)
                X=this.factor_Q10_Na_kin*90.83*exp(-(V+60)/20);
            end
        end %function
        
        function X = alpha_n(this,V)
           temperature = this.T;
            X=0;
            if (temperature<24)
                X=(this.factor_Q10_K_kin*-0.0984*(V+32.5))./(exp(-0.1*(V+32.5))-1);
            elseif(temperature>30)
                X=(this.this.factor_Q10_K_kin*-0.09575*(V+37))./(exp(-0.1*(V+37))-1);
            end
        end %function
        
        function X = beta_n(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=this.factor_Q10_K_kin*1.969*exp(-(V+58.5)/76);
            elseif(temperature>30)
                X=this.factor_Q10_K_kin*1.915*exp(-(V+47)/80);
            end
        end %function
        
        function X = alpha_c(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=(this.factor_Q10_Ca_kin*-1.4*(V+15))./(exp(-0.1*(V+15))-1);
            elseif(temperature>30)
                X=(this.factor_Q10_Ca_kin*-1.362*(V+13))./(exp(-0.1*(V+13))-1);
            end
        end %function
        
        function X = beta_c(this,V)
            temperature = this.T;
            X=0;
            if (temperature<24)
                X=this.factor_Q10_Ca_kin*46.68*exp(-(V+40)/18);
            elseif(temperature>30)
                X=this.factor_Q10_Ca_kin*45.41*exp(-(V+38)/18);
            end
        end %function
        
        function X = gKCaBar(this,Ca2i)
            X = this.gKCa*((Ca2i./this.CaDiss).^2)./(1+(Ca2i/this.CaDiss).^2);
        end % function
        
        function [T, Y] = postProcessing(this,T,Y)           
           h.figure = figure;
           plotaxis.Xmin=0; plotaxis.Xmax=this.simulationTime; plotaxis.Ymin=-90; plotaxis.Ymax=50;
           plot(T,Y(:,1,1));
           axis([plotaxis.Xmin plotaxis.Xmax plotaxis.Ymin plotaxis.Ymax]);
            
            filename=[  '/tmp/membranPotential/'  'fcn2010_cat20' '+' num2str(this.stimulationDesc(3,1)) 'ms' '+0uA'];
            filename=strrep(filename, '.', '_');
            
            if ~isequal(exist('/tmp/membranPotential', 'dir'),7) % 7 = directory.
                mkdir('/tmp/membranPotential');
            end %  if ~isequal(exist('/tmp/membranPotential', 'dir'),7)
            
            close all;
        end %function postProcessing
        
    end %method static
end %classdef
