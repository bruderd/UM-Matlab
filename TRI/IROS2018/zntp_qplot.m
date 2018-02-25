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
    patch(vx(zntp)*1e3, rad2deg(vy(zntp)), 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
%     plot(ztmax(3,:),ztmax(6,:),'r*');
    
%     quiver((zeros(1,params.num) + q(1))*1e3, rad2deg(zeros(1,params.num) + q(2)), ztmax(3,:)*1e3, rad2deg(ztmax(6,:)), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Color', 'r', 'LineWidth', 1)
    % Plot the zonotope generating vectors in separate colors
    quiver(q(1)*1e3, rad2deg(q(2)), ztmax(3,1)*1e3, rad2deg(ztmax(6,1)), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [27,158,119]/255, 'LineWidth', 2) % red [227,26,28]/255
    quiver(q(1)*1e3, rad2deg(q(2)), ztmax(3,2)*1e3, rad2deg(ztmax(6,2)), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [217,95,2]/255, 'LineWidth', 2)  % yellow [254,204,92]/255
    quiver(q(1)*1e3, rad2deg(q(2)), ztmax(3,3)*1e3, rad2deg(ztmax(6,3)), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', '.', 'Color', [117,112,179]/255, 'LineWidth', 2)  % orange [253,141,60]/255
    hold off
    
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
xlabel('Extension, $\Delta l$ (mm)', 'Interpreter', 'LaTex')
ylabel('Rotation, $\Delta \phi$ (deg)', 'Interpreter', 'LaTex')
hold off

end

