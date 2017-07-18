function vert = verticesEff( eulerAng, coordinates, params )
%Plots the end effector in orientation specified by the Euler Angles
%   Detailed explanation goes here

R = R10(eulerAng(1), eulerAng(2), eulerAng(3)); 

% defining the vertices of the end effector prism
thickness = params.thickness;   % unused for now, not actually necessary for visualization
width = params.width/2;

vert = zeros(3,8);
vert(:,1) = R*[1; 1; 0] * width + coordinates;
vert(:,2) = R*[1; -1; 0] * width + coordinates;
vert(:,3) = R*[-1; -1; 0] * width + coordinates;
vert(:,4) = R*[-1; 1; 0] * width + coordinates;
vert(:,5) = R*[width; width; thickness] + coordinates;
vert(:,6) = R*[width; -width; thickness] + coordinates;
vert(:,7) = R*[-width; -width; thickness] + coordinates;
vert(:,8) = R*[-width; width; thickness] + coordinates;

end

