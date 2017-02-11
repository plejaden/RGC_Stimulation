%% creates a vector with radius of each compartment
function compNeuronRadiusVector = createNeuronRadiusVector(geometryProperties,radiusNerve)

% display(' densityVectorFactory');
compNeuronRadiusVector= zeros(geometryProperties.numberOfCompartments,1);
for i=1:size( geometryProperties.compFeatures,1)
    thisCompartmentFeatureString =  geometryProperties.compFeatures{i,1};
    for j=1:size(geometryProperties.desc.neuronPartNames,1)
        thisCompartmentName = geometryProperties.desc.neuronPartNames(j,2);
        if (strcmp(thisCompartmentFeatureString,thisCompartmentName))
            compNeuronRadiusVector(i)= radiusNerve(j);
        end;
    end
    
end
end
