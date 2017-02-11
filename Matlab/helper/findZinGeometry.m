function [ z ] = findZinGeometry( x,y,neuronCoords,defaultZ )
z = 0;
for i=1:size(neuronCoords.name,1)
    thisY=neuronCoords.y(i);
    thisZ=neuronCoords.z(i);
    
    if y<=thisY && z==0;
        z = thisZ;
    end
end

if z==0
    z = defaultZ;
end

end

