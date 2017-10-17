function vert = verticesEff( xeul, xcart, params )
%Plots the end effector in orientation specified by the Euler Angles
%   Detailed explanation goes here

R = R10(xeul); 

% defining the vertices of the end effector prism
thickness = params.thickness;   % unused for now, not actually necessary for visualization
width = params.width/2;

vert = zeros(3,8);
vert(:,1) = R*[1; 1; 0] * width + xcart(1:3,1);
vert(:,2) = R*[1; -1; 0] * width + xcart(1:3,1);
vert(:,3) = R*[-1; -1; 0] * width + xcart(1:3,1);
vert(:,4) = R*[-1; 1; 0] * width + xcart(1:3,1);
vert(:,5) = R*[width; width; thickness] + xcart(1:3,1);
vert(:,6) = R*[width; -width; thickness] + xcart(1:3,1);
vert(:,7) = R*[-width; -width; thickness] + xcart(1:3,1);
vert(:,8) = R*[-width; width; thickness] + xcart(1:3,1);

end

