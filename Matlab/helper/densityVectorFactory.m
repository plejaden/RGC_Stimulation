%% this function matches compartment names to an vector of gBarValues
function compNeuronTypeVector = densityVectorFactory(features,compartmentNames, gBarValues)

% display(' densityVectorFactory');
compNeuronTypeVector= zeros(size(features,1),1);
for i=1:size(features,1)
    thisCompartmentFeatureString = features{i,1};
    for j=1:size(compartmentNames,2)
        thisCompartmentName = compartmentNames(1,j);
        if (strcmp(thisCompartmentFeatureString,thisCompartmentName))
            compNeuronTypeVector(i)= gBarValues(1,j);
        end;
    end
    
end
end