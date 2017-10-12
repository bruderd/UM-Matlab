function drawModule_wInputs( xeul, params, t, u )
%Visualizes the end effector position in 3D. This does not create a figure
%so is useful for animations that must iterative call a plot function.
%   Detailed explanation goes here

if nargin == 2
    t = 69;
    u = 69*ones(params.m,1);
end

L = params.Lspine;

xcart = euler2cart(xeul, params);
R = R10(xeul);

% Get the coordinates of the vertices of the end effector
vertEff = verticesEff(xeul, xcart, params);

% Define the coordinates of the vertices of the top block
width = params.width/2;
vertTop = [[width; width; -L], [width; -width; -L], [-width; -width; -L], [-width; width; -L]];
    
% Define the central spine
res = 100;
dl = L/res;
m = (R(:,3) - [0, 0, 1]') * (1/L);
b = [0; 0; 1];
spine(:,1) = [0; 0; -L];
for i = 2:res
    l = i*dl;
    spine(:,i) = spine(:,i-1) + (m*l + b) * dl;
end


% Convert pressures to KPa for display purposes
Pkpa = [u(1); u(2); u(3); u(4)] * 10^(-3);

%% Create the plot

% Draw the module geometry in 3D-space
subplot(2,1,1)
color = [189 215 231]./256;
% title(['$t = $' num2str(t,'%.2f') ':  $P_1 = $' num2str(Pkpa(1),'%.0f') ', $P_2 = $' num2str(Pkpa(2),'%.0f') ', $P_3 = $' num2str(Pkpa(3),'%.0f') ', $P_4 = $' num2str(Pkpa(4),'%.0f') ' (kPa)'], 'Interpreter', 'LaTex')
title(['$t = $' num2str(t,'%.2f') ' (s)'], 'Interpreter', 'LaTex')
hold on
patch(vertEff(1,1:4), vertEff(2,1:4), vertEff(3,1:4), color)
patch(vertEff(1,5:8), vertEff(2,5:8), vertEff(3,5:8), color)
patch([vertEff(1,1:2), vertEff(1,6), vertEff(1,5)], [vertEff(2,1:2), vertEff(2,6), vertEff(2,5)], [vertEff(3,1:2), vertEff(3,6), vertEff(3,5)], 'r')
patch([vertEff(1,2:3), vertEff(1,7), vertEff(1,6)], [vertEff(2,2:3), vertEff(2,7), vertEff(2,6)], [vertEff(3,2:3), vertEff(3,7), vertEff(3,6)], 'b')
patch([vertEff(1,3:4), vertEff(1,8), vertEff(1,7)], [vertEff(2,3:4), vertEff(2,8), vertEff(2,7)], [vertEff(3,3:4), vertEff(3,8), vertEff(3,7)], 'g')
patch([vertEff(1,1), vertEff(1,4), vertEff(1,8), vertEff(1,5)], [vertEff(2,1), vertEff(2,4), vertEff(2,8), vertEff(2,5)], [vertEff(3,1), vertEff(3,4), vertEff(3,8), vertEff(3,5)], 'm')
patch(vertTop(1,:), vertTop(2,:), vertTop(3,:), color)
plot3(spine(1,:), spine(2,:), spine(3,:), 'LineWidth',5)
hold off
set(gca,'zdir','reverse')
view(3)
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal
xlim([-L, L])
ylim([-L, L])
zlim([-L, 5e-2])
box on

% Plot the input pressure to each FREE
subplot(2,1,2)
bar(Pkpa)
title(['$P_1 = $' num2str(Pkpa(1),'%.0f') ', $P_2 = $' num2str(Pkpa(2),'%.0f') ', $P_3 = $' num2str(Pkpa(3),'%.0f') ', $P_4 = $' num2str(Pkpa(4),'%.0f') ' (kPa)'], 'Interpreter', 'LaTex')
ylim([0, 3.1e6*1e-3])
ylabel('Input Pressure (kPa)')


end
