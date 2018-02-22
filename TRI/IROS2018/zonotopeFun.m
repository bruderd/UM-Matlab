function [zntp, vx, vy] = zonotopeFun(x, params)

% Set the displacement of the end effector
% x = [0, 0, 0, 0, 0, 0]';
q = x2q(x, params);

% Calculate the maximum Z = [F,M]' for each FREE
Zmax = maxZ(q, params);

% Convert Zmax to end effector coordinates for each FREE
zetamax = Z2zeta_i(Zmax, params);

% Generate the zonotope
ztmax = zeros(6,params.num);
for i = 1:params.num
    ztmax(:,i) = zetamax(6*(i-1)+1: 6*i, 1);   % stack them horizontially so genZonotope can read them
end
[zntp, vx, vy] = genZonotope(ztmax);

% Plot the force zonotope
figure
hold on
% quiver(zeros(params.num,1), zeros(params.num,1), ztmax(3,:)',ztmax(6,:)')
quiver(zeros(1,params.num), zeros(1,params.num), ztmax(3,:),ztmax(6,:), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', 'o', 'Color', 'r', 'LineWidth', 2)
patch(vx(zntp), vy(zntp), 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
plot(ztmax(3,:),ztmax(6,:),'r*');
% plot(vx(zntp), vy(zntp), 'b-');
xlim([-15 15]);
ylim([-0.1 0.1]);
xL = xlim;
yL = ylim;
line([0 0], yL, 'color', 'k');  %x-axis
line(xL, [0 0], 'color', 'k');  %y-axis
grid on
box on
% Commented out, will just label in latex...
% xlabel('Force (N)')
% ylabel('Moment (N-m)')
hold off

