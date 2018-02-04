% zonotopeMain.m

% Set the values of all of the parameters
Gama = deg2rad([20, -20, 85]);
R = 0.01 * ones(1,3);
L = 0.1 * ones(1,3);
d = zeros(3,3);
p = [0,0,0 ; 0,0,0 ; 1,1,1];
Pmax = [100, 100, 1000];
params = setParams(Gama, R, L, d, p, Pmax);

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

