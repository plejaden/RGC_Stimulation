%% creates a vector with information, which neuronal part a specific compartmen belongs to.
function compNeuronTypeVector = createCompNeuronTypeVector(geometryProperties,measurementPoints,neuronCoords)
%display(' createCompNeuronTypeVektor');

compNeuronTypeVector=cell(length(measurementPoints),1);
geometryDescription = geometryProperties.desc;

%% loop over all measurment-points
for i=1:length(measurementPoints)
    haystackX=measurementPoints(i,1);
    haystackY=measurementPoints(i,2);
    haystackZ=measurementPoints(i,3);
    
    %% loop over all neural parts
    for j=1:size(neuronCoords.x,1)
        needleX=neuronCoords.x(j);
        needleY=neuronCoords.y(j);
        needleZ=neuronCoords.z(j);
        neuronPartName= neuronCoords.name(j);
        
        if (needleY<=haystackY )
            neuronPartName = descriptionMatch(geometryDescription,neuronPartName);
            compNeuronTypeVector{i}=neuronPartName;
        end
        
        
    end %  for j=1:size(neuronCoords.x,1)-1
    
end %for i=1:length(measurementPoints)
%display(compNeuronTypeVector);
end


function name = descriptionMatch(geometryDescription,featureName)
%display(' descriptionMatch');
name='';
thisNeuronPartNames=geometryDescription.neuronPartNames;
%% loop over all neural parts
for j=1:size(thisNeuronPartNames,1)
    
    if (strcmp(thisNeuronPartNames(j,1),cell2mat(featureName{1})))
        name= thisNeuronPartNames(j,2);
        return;
    end
    
    
end %  for j=1:size(neuronCoords.x,1)-1
end
