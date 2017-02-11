function [ desc ] = description_saiful_SomaIsAxon_epiretinal_075um_s_development( )
desc.neuronPartNames={'sph1' 'Soma';'cyl11' 'IS';'cyl9' 'Axon';'cyl10' 'SOCB'};
desc.electrodeGeometryNames={'cyl2';'cyl3';'cyl4';'cyl5';'cyl6';'cyl7';'cyl8'};
desc.variables={'normalDistA=0.075/2';'normalDistC=normalDistA/sin(60*2*pi/360)';'normalDistB=normalDistC*0.5';'distancePower=10^-3';'electrodeHeight=10*10^-6';'electrodeZpos=315*10^-6';'electrodeEdge=20*10^-6';'neuronYpos=1.98*10^-3'};
end
