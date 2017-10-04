function plotPos_array( P1, P2, P3, P4, params)
%Creates an array of plots of the position of the end effector for
%different pressure combinations
%   Detailed explanation goes here

% Convert pressures to KPa for display purposes
P1kpa = P1 * 10^(-3);
P2kpa = P2 * 10^(-3);
P3kpa = P3 * 10^(-3);
P4kpa = P4 * 10^(-3);

fig1 = plotPos(P1, params);
fig2 = plotPos(P2, params);
fig3 = plotPos(P3, params);
fig4 = plotPos(P4, params);

figure
h(1) = subplot(2,2,1);
title(['$P_1 = $' num2str(P1kpa(1)) ', $P_2 = $' num2str(P1kpa(2)) ', $P_3 = $' num2str(P1kpa(3)) ', $P_4 = $' num2str(P1kpa(4)) ' (kPa)'], 'Interpreter', 'LaTex')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(3)
set(gca,'zdir','reverse')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-params.Lspine, inf])
box on

h(2) = subplot(2,2,2);
title(['$P_1 = $' num2str(P2kpa(1)) ', $P_2 = $' num2str(P2kpa(2)) ', $P_3 = $' num2str(P2kpa(3)) ', $P_4 = $' num2str(P2kpa(4)) ' (kPa)'], 'Interpreter', 'LaTex')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(3)
set(gca,'zdir','reverse')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-params.Lspine, inf])
box on

h(3) = subplot(2,2,3);
title(['$P_1 = $' num2str(P3kpa(1)) ', $P_2 = $' num2str(P3kpa(2)) ', $P_3 = $' num2str(P3kpa(3)) ', $P_4 = $' num2str(P3kpa(4)) ' (kPa)'], 'Interpreter', 'LaTex')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(3)
set(gca,'zdir','reverse')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-params.Lspine, inf])
box on

h(4) = subplot(2,2,4);
title(['$P_1 = $' num2str(P4kpa(1)) ', $P_2 = $' num2str(P4kpa(2)) ', $P_3 = $' num2str(P4kpa(3)) ', $P_4 = $' num2str(P4kpa(4)) ' (kPa)'], 'Interpreter', 'LaTex')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(3)
set(gca,'zdir','reverse')
axis equal
xlim([-10e-2, 10e-2])
ylim([-10e-2, 10e-2])
zlim([-params.Lspine, inf])
box on

copyobj(allchild(get(fig1,'CurrentAxes')),h(1));
copyobj(allchild(get(fig2,'CurrentAxes')),h(2));
copyobj(allchild(get(fig3,'CurrentAxes')),h(3));
copyobj(allchild(get(fig4,'CurrentAxes')),h(4));

end

