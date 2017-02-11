classdef fm1997_singleComp < basicModel
% Last edit: 2016-03-10 16:36:45 
    
    properties
        simulationTime=10;
        t0; timeout=300; %s
        
        gasConstant=8.31441; % (J/mol*K)
        Faraday=96485; %c/mol
        RTdivF; %mV
        Temp=293.15; % 20 Degree Celsius
        area;
        volume;
        gNa=50.0;%(mS/cm^2)
        gCa=2.2;% (mS/cm2)
        gK=12;% (mS/cm2)
        gKA=36;% (mS/cm2)
        gKCa=0.05;% (mS/cm2)
        gL=0.05; %(mS/cm2)
        C = 1; %(mF/cm2)
        Vrest=-65; % (mV)
        VNa = +35; %(mV)
        VK = -75; % (mV)
        VL=-60; %(mV) Vleakage = -60 to -65 mV
        Rn = 10^9;%(ohm) = 1 GOhm
        
        gA=36; %
        r;compartmentLength;
        CaRes=10^-7*10^3; % in M/dm3 => mM/cm3
        CaDiss=10^-6*10^3; % in M/dm3 => mM/cm3
        TauCA=1.5; %
        Cao=1.8; % mM
        
        %% simulation
        iPause=0.5;
        iStim=0; %(uA/cm^2) Amount of current stimulation for 0.12msec
        iStimTime=0.120;
        
        % biphasic stimulation: for this time, the mentioned multiplyer will be applied %
        % format: (time in [ms], pot multiplicator [1] ; )
        % [ puls1,pot1 ; puls2,pot2  ; puls3,pot3 ; puls4,pot4]
        % exmaple: [0.5,0.15,0.2,0.1; 0,-1,0,+1]
        
        % intracellular stimulation:
        stimulationDescInt=[0.5,0; 0.2,1; 4.5,0; 0.2,1];
        
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
        function this = fm1997_singleComp(geometryProperties)
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
            
            l=geometryProperties.numberOfCompartments; % get the number of compartments from the voltage array%
            
            %% current analysis
            currentObserver.compartment=1; currentObserver.time = 0; currentObserver.INa = 0; currentObserver.IK = 0;
            currentObserver.IKA = 0; currentObserver.IKCa = 0; currentObserver.ICa = 0; currentObserver.IL = 0; currentObserver.stimulusCurrent= 0;
            
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
            thisINa=zeros(compL,1); thisIK=zeros(compL,1); thisIKCa=zeros(compL,1);
            thisICa=zeros(compL,1); thisIL=zeros(compL,1); thisIKA=zeros(compL,1);
            
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
            %             plot(T,Y(:,1,1));
            %             xlabel('t [ms]'); ylabel('membran potential [mV]');
            filename=[  '/tmp/membranPotential/'  'fm1997_singlecomp' '+' num2str(this.stimulationDesc(3,1)) 'ms'];
            filename=strrep(filename, '.', '_');
            
            if ~isequal(exist('/tmp/membranPotential', 'dir'),7) % 7 = directory.
                mkdir('/tmp/membranPotential');
            end
            
            save (filename);
            %             mkdir('/tmp/membranPotential');
            %             myStyle = hgexport('factorystyle'); myStyle.Format = 'png'; myStyle.Width = 3; myStyle.Height = 1; myStyle.Resolution = 300; myStyle.Units = 'inch'; myStyle.FixedFontSize = 12;
            %             hgexport(gcf, [ '/tmp' '/' 'membranPotential' '/' filename '.png' ], myStyle, 'Format', 'png');
            %             close all;
        end %function postProcessing
        
    end %method static
    
end %classdef
