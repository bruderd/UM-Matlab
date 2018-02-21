function zntp_qplot( qpoints, scale, params)
%zntp_qplot: Plots the force zonotope for different values of s,w on same
%plot
%   xpoints: Points at which to calculate the zonotope, [s1, s2, ...; w1, w2...]
%   params: system parameters
%   testParams: testing parameters, like maximum pressure

xscale = scale(1);  % scaling between forces and displacements for visual clarity
yscale = scale(2);  % scaling between forces and displacements for visual clarity

figure
for i = 1: length(qpoints(1,:))
    q = qpoints(:,i);
    [zntp, ztmax, vx, vy] = zonotopeFun_noPlot(q, params);
    vx = xscale*vx + q(1);
    vy = yscale*vy + q(2);
    ztmax(3,:) = ztmax(3,:) * xscale;
    ztmax(6,:) = ztmax(6,:) * yscale;
   
    % Plot the force zonotope
    hold on
%     plot(q(1), q(2), 'r*')
    quiver(zeros(1,params.num) + q(1), zeros(1,params.num) + q(2), ztmax(3,:),ztmax(6,:), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Color', 'r', 'LineWidth', 1)
    patch(vx(zntp), vy(zntp), 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
%     plot(ztmax(3,:),ztmax(6,:),'r*');
end

% change/specify plot features
hold on
% xlim([-15 15]);
% ylim([-0.1 0.1]);
xL = xlim;
yL = ylim;
line([0 0], yL, 'color', 'k');  %x-axis
line(xL, [0 0], 'color', 'k');  %y-axis
grid on
box on
xlabel('Extension, $s$ (m)', 'Interpreter', 'LaTex')
ylabel('Rotation, $w$ (rad)', 'Interpreter', 'LaTex')
hold off

end

