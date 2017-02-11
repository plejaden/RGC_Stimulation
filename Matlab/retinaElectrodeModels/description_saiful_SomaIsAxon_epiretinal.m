function [ desc ] = description_saiful_SomaIsAxon_epiretinal( )
desc.neuronPartNames={'sph1' 'Soma';'cyl2' 'IS';'cyl1' 'Axon'};
desc.electrodeGeometryNames={'blk4';'blk5';'blk6';'blk7';'blk8';'blk9';'blk10'};
desc.variables={'normalDistA=0.5';'normalDistC=normalDistA/sin(60*2*pi/360)';'normalDistB=normalDistC*0.5';'distancePower=10^-3';'electrodeHeight=10*10^-6';'electrodeZpos=390*10^-6'};
end
