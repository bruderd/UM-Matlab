function fig = plotPos( P, params )
%Visualizes the end effector position in 3D
%   Detailed explanation goes here

L = params.Lspine;

eulerAng = findEulerAng(P, params);
coordinates = orient2pos(eulerAng, params);
R = R10(eulerAng(1), eulerAng(2), eulerAng(3));

% Get the coordinates of the vertices of the end effector
vertEff = verticesEff(eulerAng, coordinates, params);

% Define the coordinates of the vertices of the top block
width = params.width/2;
vertTop = [[width; width; -L], [width; -width; -L], [-width; -width; -L], [-width; width; -L]];
    
% Define the central spine
res = 100;
dl = L/res;
m = R*[0; 0; 1/L] - [0; 0; 1/L];
b = [0; 0; 1];
spine(:,1) = [0; 0; -L];
for i = 2:100
    l = i*dl;
    spine(:,i) = spine(:,i-1) + (m*l + b) * dl;
end


% Create the plot
color = [189 215 231]./256;

fig = figure('Name','End Effector Position');
hold on
patch(vertEff(1,1:4), vertEff(2,1:4), vertEff(3,1:4), color)
patch(vertEff(1,5:8), vertEff(2,5:8), vertEff(3,5:8), color)
patch([vertEff(1,1:2), vertEff(1,6), vertEff(1,5)], [vertEff(2,1:2), vertEff(2,6), vertEff(2,5)], [vertEff(3,1:2), vertEff(3,6), vertEff(3,5)], color)
patch([vertEff(1,2:3), vertEff(1,7), vertEff(1,6)], [vertEff(2,2:3), vertEff(2,7), vertEff(2,6)], [vertEff(3,2:3), vertEff(3,7), vertEff(3,6)], color)
patch([vertEff(1,3:4), vertEff(1,8), vertEff(1,7)], [vertEff(2,3:4), vertEff(2,8), vertEff(2,7)], [vertEff(3,3:4), vertEff(3,8), vertEff(3,7)], color)
patch([vertEff(1,1), vertEff(1,4), vertEff(1,8), vertEff(1,5)], [vertEff(2,1), vertEff(2,4), vertEff(2,8), vertEff(2,5)], [vertEff(3,1), vertEff(3,4), vertEff(3,8), vertEff(3,5)], color)
patch(vertTop(1,:), vertTop(2,:), vertTop(3,:), color)
plot3(spine(1,:), spine(2,:), spine(3,:), 'LineWidth',5)
hold off
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-params.Lspine, inf])
box on


end

