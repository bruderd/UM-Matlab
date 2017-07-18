%plots_KvP_thetavP.m
% Plot values of theta and K at a bunch of different combinations of
% pressures
%
% Must run setParams.m before this script

% set the pressure step sizes based on the min/max from setParams_2parallel.m
dP1 = (params.Prange1(2) - params.Prange1(1)) / params.steps1;
dP2 = (params.Prange2(2) - params.Prange2(1)) / params.steps2;

% initializations
P1 = zeros(1,params.steps1);
P2 = zeros(1,params.steps2);
thetaSS = zeros(params.steps1, params.steps2);
K = zeros(params.steps1, params.steps2);

for i = 1:params.steps1
    
    P1(i) = i * dP1;
    
    for j = 1:params.steps2
        
        P2(j) = j * dP2;
        
        u = [P1(i); P2(j)];
        
        thetaSS(i,j) = find_thetaSS(u, params);
        
        K(i,j) = find_K(thetaSS(i,j), u, params);

    end
    
end

%% Generate plots for theta vs. P and K vs. P

% change units to make plots easier to understand
P_kPa = [P1; P2] * 1e-3;

[X,Y] = meshgrid(P_kPa(1,:), P_kPa(2,:));

figure
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)

% plot P vs. thetaSS (colormap)
subplot(1,2,1)
surf(X,Y,thetaSS)
view(2)
sbar = colorbar;
ylabel(sbar, '$\theta_{ss}$ (rad)', 'fontsize', 20, 'Interpreter', 'LaTex')
title('Rotation', 'fontsize', 20)
xlabel('$P_{1}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
ylabel('$P_{2}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
zlabel('s (mm)', 'fontsize', 20, 'Interpreter', 'LaTex')

% plot P vs. K (colormap)
subplot(1,2,2)
surf(X,Y,K)
view(2)
wbar = colorbar;
ylabel(wbar, '$K$ (N/m)', 'fontsize', 20, 'Interpreter', 'LaTex')
title('Stiffness', 'fontsize', 20)
xlabel('$P_{1}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
ylabel('$P_{2}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
zlabel('$\omega$ (degrees)', 'fontsize', 20, 'Interpreter', 'LaTex')

