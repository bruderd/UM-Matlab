function fig = plotOrient( P, params )
%Visualizes the end effector orientation in 3D
%   Detailed explanation goes here

% Rename relevant parameters
a = [params.attach1(1), params.attach2(1), params.attach3(1), params.attach4(1)];
b = [params.attach1(2), params.attach2(2), params.attach3(2), params.attach4(2)];

% Determine the Euler angles and corresponding rotation matrix
eulerAng = findEulerAng(P, params);
R = R10(eulerAng(1), eulerAng(2), eulerAng(3)); 

% defining the vertices of the end effector prism
thickness = params.thickness;
width = params.width/2;
top = [R*[1; 1; 0], R*[1; -1; 0], R*[-1; -1; 0], R*[-1; 1; 0]] * width;
bottom = [R*[1; 1; thickness], R*[1; -1; thickness], R*[-1; -1; thickness], R*[-1; 1; thickness]] * width;
side1 = [R*[1; 1; 0], R*[1; 1; thickness], R*[1; -1; thickness], R*[1; -1; 0]] * width;
side2 = [R*[1; 1; 0], R*[1; 1; thickness], R*[-1; 1; thickness], R*[-1; 1; 0]] * width;
side3 = [R*[-1; -1; 0], R*[-1; -1; thickness], R*[1; -1; thickness], R*[1; -1; 0]] * width;
side4 = [R*[-1; -1; 0], R*[-1; -1; thickness], R*[-1; 1; thickness], R*[-1; 1; 0]] * width;

color = [189 215 231]./256;

fig = figure('Name','End Effector Orientation');
patch(top(1,:), top(2,:), top(3,:), color)
patch(bottom(1,:), bottom(2,:), bottom(3,:), color)
patch(side1(1,:), side1(2,:), side1(3,:), color)
patch(side2(1,:), side2(2,:), side2(3,:), color)
patch(side3(1,:), side3(2,:), side3(3,:), color)
patch(side4(1,:), side4(2,:), side4(3,:), color)
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-inf, inf])


end

