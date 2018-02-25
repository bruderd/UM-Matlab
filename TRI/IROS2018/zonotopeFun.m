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
% % Plot all zonotope generating vectors in one color
% quiver(zeros(1,params.num), zeros(1,params.num), ztmax(3,:),ztmax(6,:), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', 'o', 'Color', 'r', 'LineWidth', 2)

% Fill in the zonotope
patch(vx(zntp), vy(zntp), 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
% plot(ztmax(3,:),ztmax(6,:),'r*');         % ends of the zonotope generating vectors
% plot(vx(zntp), vy(zntp), 'k-', 'LineWidth', 1);     % zonotope outline

% Plot the zonotope generating vectors in separate colors
quiver(0, 0, ztmax(3,1),ztmax(6,1), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [27,158,119]/255, 'LineWidth', 5) % red [227,26,28]/255
quiver(0, 0, ztmax(3,2),ztmax(6,2), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [217,95,2]/255, 'LineWidth', 5)  % yellow [254,204,92]/255
quiver(0, 0, ztmax(3,3),ztmax(6,3), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [117,112,179]/255, 'LineWidth', 5)  % orange [253,141,60]/255

xlim([-15 15]);
ylim([-0.1 0.1]);
xL = xlim;
yL = ylim;
line([0 0], yL, 'color', 'k');  %x-axis
line(xL, [0 0], 'color', 'k');  %y-axis
grid on
box on
% Commented out, will just label in latex...
xlabel('Force, $F^{\hat{z}_e}$ (N)', 'Interpreter', 'LaTex')
ylabel('Moment, $M^{\hat{z}_e}$ (N-m)', 'Interpreter', 'LaTex')
hold off

