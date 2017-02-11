%% 2015-04-23 17:18:19
classdef fm1997_multiComp < basicModel
    
    properties
        simulationTime=10;
        t0; timeout=300; %s
        
        gasConstant=8.31441; % (J/mol*K)
        Faraday=96485; %c/mol
        RTdivF; %mV
        Temp=293.15; % 20 Degree Celsius
        area;
        volume;
        gNa; gK; gCa; gKA; gKCa; gL;
        
        %% channel densities soma (Equivalent Cylinder EC 1um)
        gNaSoma=50;%(mS/cm^2)
        gKSoma=12;% (mS/cm2)
        gCaSoma=1;% (mS/cm2)
        gKASoma=36; % (mS/cm2)
        gKCaSoma=0.065;% (mS/cm2)
        gLSoma=0.003; % (mS/cm2)
        
        
        %% channel densities axon (Equivalent Cylinder EC 1um)
        gNaAxon=70;%(mS/cm^2)
        gKAxon=70;% (mS/cm2)
        gCaAxon=04;% (mS/cm2)
        gKAAxon=0; % (mS/cm2)
        gKCaAxon=0.065;% (mS/cm2)
        gLAxon=0.003; % (mS/cm2)
        
        
        %% elibrium potentials
        VNa = 35; %(mV)
        VK = -101.34; % (mV)
        VL= -30; %(mV)
        
        %% Calcium properties
        CaRes=10^-7*10^3; % in M/dm3 => mM/cm3
        CaDiss=10^-6*10^3; % in M/dm3 => mM/cm3
        Cao=1.8; % mM
        
        %% misc properties
        Vrest= -70; % (mV)
        Rn = 150*10^6; %(ohm) = 150 MOhm, input resistance
        gA=5.4; %
        r;l;
        C = 1; %(mF/cm2)
        TauCA=1.5; %
        
        %% simulation
        iPause=0.5;
        iStim; %(uA/cm^2) Amount of current stimulation for 0.12msec
        iStimTime=0.120;
        
        % biphasic stimulation: for this time, the mentioned multiplyer will be applied %
        % format: (time in [ms], pot multiplicator [1] ; )
        % [ puls1,pot1 ; puls2,pot2  ; puls3,pot3 ; puls4,pot4]
        % exmaple: [0.5,0.15,0.2,0.1; 0,-1,0,+1]
        
        % intracellular stimulation:
        stimulationDescInt=[0.5,0; 0.2,1; 4,0; 0.2,1];
        
        % extracellular stimulation:
        stimulationDescExt=[0.5,0; 0.2,1; 1,0; 0.4,-0.5]; % with electrode uncharge
        % stimulationDescExt=[0.5,0; 0.2,1; 1,0; 0.4,0]; % without electrode uncharge
        
        startParameters;
        numberOfStateVariables=7; % h m n c a hA Cai
        Ve=0;
        rhoFiber = 1.1; %[Ohm*m]
        stimulationDesc; ohmAx;
        
    end %properties
    
    methods
        
        %% The standard constructor sets the temperature of the model.
        function this = fm1997_multiComp(geometryProperties)
            global currentObserver;
            
            %prevents error if t0 is not initialized
            if (isempty(this.t0)) this.t0=clock; end
            
            ms = round(etime(clock,this.t0) * 1000);
            if(ms>this.timeout*1000)
                return;
            end
            
            this.Ve=geometryProperties.externalVoltages;
            this.iStim=geometryProperties.iStim;
            this.r=geometryProperties.compartmentRadius*100; % m => cm
            this.compartmentLength=geometryProperties.distanceBetrweenCompartments*100; % m => cm
            this.area=this.r.^2*pi;   % area in cm2
            this.volume=geometryProperties.distanceBetrweenCompartments*this.area;
            
            if (this.iStim==0)
                display('using extracellular stimulation function');
                this.stimulationDesc=this.stimulationDescExt;
            else
                display('using intercellular stimulation function');
                this.stimulationDesc=this.stimulationDescInt;
            end % (this.iStim~=0)
            
            % channel densities vector
            this.gNa = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gNaSoma,this.gNaSoma,this.gNaAxon]);
            this.gK = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gKSoma,this.gKSoma,this.gKAxon]);
            this.gCa = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gCaSoma,this.gCaSoma,this.gCaAxon]);
            this.gKA = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gKASoma,this.gKASoma,this.gKAAxon]);
            this.gKCa = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gKCaSoma,this.gKCaSoma,this.gKCaAxon]);
            this.gL = densityVectorFactory(geometryProperties.compFeatures,{'Soma' 'IS' 'Axon'},[this.gLSoma,this.gLSoma,this.gLAxon]);
            
            %% current analysis
            currentObserver.compartment=1; currentObserver.time = 0; currentObserver.INa = 0; currentObserver.IK = 0;
            currentObserver.IKA = 0; currentObserver.IKCa = 0; currentObserver.ICa = 0; currentObserver.IL = 0; currentObserver.stimulusCurrent= 0;
            
            l=geometryProperties.numberOfCompartments; % get the number of compartments from the voltage array%
            V0=repmat(this.Vrest,1,l)';
            this.startParameters=[V0 this.h0(l) this.m0(l) this.n0(l) this.c0(l) this.A0(l) this.hA0(l) this.Cai0(l)];
            this.RTdivF=(this.gasConstant*this.Temp)/this.Faraday;
            
        end
        
        
        %% the differential ewuations
        function X = differentialEq(this,time,x)
            
            global currentObserver;
            
            % slice up the column vector to a 'compL'x8 array
            compL=length(x)/(this.numberOfStateVariables+1);
            V=x(1:compL); h=x(compL+1:2*compL); m=x(2*compL+1:3*compL); n=x(3*compL+1:4*compL);
            c=x(4*compL+1:5*compL); a=x(5*compL+1:6*compL); hA=x(6*compL+1:7*compL); Cai=abs(x(7*compL+1:8*compL));
            thisVe=this.Ve;
            
            %% get the stimulus current for each compartment, if the nerve gets stimulated on the Soma%
            stimulusCurrent=this.linkCompartments(this,compL,time,V,this.stimulationDesc,this.iStim, thisVe); % current in mA %
            
            %% allocate memory
            X=zeros(compL,8);
            
            %% Main equation
            thisINa=this.INa(V,m,h);
            thisIK=this.IK(V,n);
            thisIKA=this.IKA(V,hA,a);
            thisIKCa=this.IKCa(V,Cai);
            thisICa=this.ICa(V,c,Cai);
            thisIL=this.IL(V);
            
            currentObserver.time = [currentObserver.time; time]; currentObserver.INa = [currentObserver.INa; thisINa(currentObserver.compartment)];
            currentObserver.IK = [currentObserver.IK; thisIK(currentObserver.compartment)];  currentObserver.IKA = [currentObserver.IKA; thisIKA(currentObserver.compartment)];
            currentObserver.IKCa = [currentObserver.IKCa; thisIKCa(currentObserver.compartment)]; currentObserver.ICa = [currentObserver.ICa; thisICa(currentObserver.compartment)];
            currentObserver.IL = [currentObserver.IL; thisIL(currentObserver.compartment)]; currentObserver.stimulusCurrent = [currentObserver.stimulusCurrent; stimulusCurrent(currentObserver.compartment)];
            
            X(:,1)=(-thisINa-thisIK-thisIKA-thisIKCa-thisICa-thisIL+stimulusCurrent)./this.C; %mV;
            
            X(:,2)= -(this.alpha_h(this.E(V))+this.beta_h(this.E(V))).*h+this.alpha_h(this.E(V)); %h
            X(:,3)= -(this.alpha_m(this.E(V))+this.beta_m(this.E(V))).*m+this.alpha_m(this.E(V)); %m
            X(:,4)= -(this.alpha_n(this.E(V))+this.beta_n(this.E(V))).*n+this.alpha_n(this.E(V)); %n
            X(:,5)= -(this.alpha_c(this.E(V))+this.beta_c(this.E(V))).*c+this.alpha_c(this.E(V)); %c
            X(:,6)= -(this.alpha_A(this.E(V))+this.beta_A(this.E(V))).*a+this.alpha_A(this.E(V)); %a
            X(:,7)= -(this.alpha_hA(this.E(V))+this.beta_hA(this.E(V))).*hA+this.alpha_hA(this.E(V)); %hA
            X(:,8)=(this.surfaceToVolumeRatio(this.r,this.compartmentLength)*(-this.ICa(V,c,Cai)./(this.Faraday*2)) * 10^-3)-((Cai-this.CaRes)./this.TauCA);%Cai
            
            %create a column vector
            X=X(:);
        end
        
        %% Helper functions
        function X = E(this,V)
            X=V;
        end %function E
        
        function X=INa(this,V,m,h) % check
            X=this.gNa.*m.^3.*h.*(V-this.VNa);
        end %function INa
        
        function X=IK(this,V, n)
            X=this.gK.*(n.^4).*(V-this.VK);
        end %function IK
        
        function X=IKA(this,V,ha,a)
            X=this.gA.*(a.^3).*ha.*(V-this.VK);
        end %function IKA
        
        function X=IKCa(this,V,Cai)
            gKCaBar = this.gKCa.*((Cai./this.CaDiss).^2)./(1+(Cai./this.CaDiss).^2);
            X=gKCaBar .* (V-this.VK);
        end %function IKCa
        
        function X=ICa(this,V,c,Cai)
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
            X=this.alpha_h(this.Vrest)./(this.alpha_h(this.Vrest)+this.beta_h(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = m0(this,l)
            X=this.alpha_m(this.Vrest)./(this.alpha_m(this.Vrest)+this.beta_m(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = n0(this,l)
            X=this.alpha_n(this.Vrest)./(this.alpha_n(this.Vrest)+this.beta_n(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = c0(this,l)
            X=this.alpha_c(this.Vrest)./(this.alpha_c(this.Vrest)+this.beta_c(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = A0(this,l)
            X=this.alpha_A(this.Vrest)./(this.alpha_A(this.Vrest)+this.beta_A(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = hA0(this,l)
            X=this.alpha_hA(this.Vrest)./(this.alpha_hA(this.Vrest)+this.beta_hA(this.Vrest));
            X=repmat(X,1,l)';
        end %function
        
        function X = Cai0(~,l)
            X=1.0000e-04; % mM
            X=repmat(X,1,l)';
        end %function
        
    end %methods
    
    methods (Static)
        function X = alpha_h(E)
            X=0.4.*exp(-(50+E)/20);
        end %function
        
        function X = beta_h(E)
            X=6./(exp(-0.1*(20+E))+1);
        end %function
        
        function X = alpha_m(E)
            X=(-0.6*(E+30))./(exp(-0.1*(E+30))-1);
        end %function
        
        function X = beta_m(E)
            X=20*exp(-(E+55)/18);
        end %function
        
        function X = alpha_n(E)
            X=(-0.02*(E+40))./(exp(-0.1*(E+40))-1);
        end %function
        
        function X = beta_n(E)
            X=0.4*exp(-(E+50)/80);
        end %function
        
        %% not inherited functions
        function X = alpha_c(E)
            X=(-0.3*(E+13))./(exp(-0.1*(E+13))-1);
        end %function
        
        function X = beta_c(E)
            X=10*exp(-(E+38)/18);
        end %function
        
        function X = alpha_A(E)
            X=(-0.006*(E+90))./(exp(-0.1*(E+90))-1);
        end %function
        
        function X = beta_A(E)
            X=0.1*exp(-(E+30)/10);
        end %function
        
        function X = alpha_hA(E)
            X=0.04*exp(-(E+70)/20);
        end %function
        
        function X = beta_hA(E)
            X=0.6./(exp(-0.1*(E+40))+1);
        end %function
        
        function [T, Y] = postProcessing(this,T,Y)
            
            h.figure = figure;
            plotaxis.Xmin=0; plotaxis.Xmax=this.simulationTime; plotaxis.Ymin=-90; plotaxis.Ymax=50;
            plot(T,Y(:,1,1));
            axis([plotaxis.Xmin plotaxis.Xmax plotaxis.Ymin plotaxis.Ymax]);
            
            filename=[ 'fm1997_multiComp' '+' num2str(this.stimulationDesc(3,1)) 'ms'];
            filename=strrep(filename, '.', '_');
            
            if ~isequal(exist('/tmp/membranPotential', 'dir'),7) % 7 = directory.
                mkdir('/tmp/membranPotential');
            end
            
            %  save (filename);
            
            %  cleanfigure;
            %  matlab2tikz([ '/tmp' '/' 'membranPotential' '/' filename '.tikz' ] , 'showInfo', false, ...
            %  'parseStrings',false, ...
            %  'standalone', false, ...
            %  'height', '\figureheight', ...
            %  'width','\figurewidth');
            %
            close all;
        end %function postProcessing
        
    end %method static
    
end %classdef
