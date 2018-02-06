function [zntp, vx, vy] = zonotopeFun(params)

% Set the displacement of the end effector
x = [0, 0, 0, 0, 0, 0]';
q = x2q(x, params);

% Calculate the maximum Z = [F,M]' for each FREE
Zmax = maxZ(q, params);

% Convert Zmax to end effector coordinates
zetamax = Z2zeta(Zmax, params);

% Generate the zonotope
ztmax = zeros(6,params.num);
for i = 1:params.num
    ztmax(:,i) = zetamax(6*(i-1)+1: 6*i, 1);   % stack them horizontially so genZonotope can read them
end
[zntp, vx, vy] = genZonotope(ztmax);

% Plot the results
figure
hold on
% quiver(zeros(params.num,1), zeros(params.num,1), ztmax(3,:)',ztmax(6,:)')
quiver(zeros(1,params.num), zeros(1,params.num), ztmax(3,:),ztmax(6,:), 'ShowArrowHead', 'off', 'AutoScaleFactor', 1, 'Marker', 'o', 'Color', 'r')
plot(ztmax(3,:),ztmax(6,:),'r*');
plot(vx(zntp), vy(zntp), 'b-');
hold off

