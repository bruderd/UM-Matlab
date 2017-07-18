% plots_qvP.m
%   Solves for q = [s, w]' over a range of pressures P = [P_l, P_r]'
%   Plots results in 3d plot and color plot

% set the pressure step sizes based on the min/max from setParams_2parallel.m
dP_l = (params.Prange_l(2) - params.Prange_l(1)) / params.steps_l;
dP_r = (params.Prange_r(2) - params.Prange_r(1)) / params.steps_r;

% initializations
P_l = zeros(1,params.steps_l);
P_r = zeros(1,params.steps_r);
s = zeros(params.steps_l, params.steps_r);
w = zeros(params.steps_l, params.steps_r);
remainder_F = zeros(params.steps_l, params.steps_r);
remainder_M = zeros(params.steps_l, params.steps_r);

for i = 1:params.steps_l
    
    P_l(i) = i * dP_l;
    
    for j = 1:params.steps_r
        
        P_r(j) = j * dP_r;
             
        q = lsqnonlin(@(q) netF(q, [P_l(i); P_r(j)], params), [0;0]);
        s(i,j) = q(1);
        w(i,j) = q(2);
        
        % calculate the remainder of the optimization
        remainder = netF(q, [P_l(i); P_r(j)], params);
        remainder_F(i,j) = remainder(1);
        remainder_M(i,j) = remainder(2);
%         P(i,j) = P_r(j);
    end
    
end


%% Generate plots for P vs. s and P vs. w

% change units to make plots easier to understand
P_kPa = [P_l; P_r] * 1e-3;
s_mm = s * 1e3;
w_deg = rad2deg(w);


figure
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)
[X,Y] = meshgrid(P_kPa(1,:), P_kPa(2,:));

% plot P vs. s (3d plot)
subplot(1,2,1)
surf(X,Y,s_mm')
% view(2)
% sbar = colorbar;
% ylabel(sbar, 's (mm)', 'fontsize', 20, 'Interpreter', 'LaTex')
title('Linear Displacement', 'fontsize', 20)
xlabel('$P_{\textit{l}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
ylabel('$P_{\textrm{r}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
zlabel('s (mm)', 'fontsize', 20, 'Interpreter', 'LaTex')

% plot P vs. w (3d plot)
subplot(1,2,2)
surf(X,Y,w_deg')
% view(2)
% wbar = colorbar;
% ylabel(wbar, '$\omega$ (degrees)', 'fontsize', 20, 'Interpreter', 'LaTex')
title('Rotational Displacement', 'fontsize', 20)
xlabel('$P_{\textit{l}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
ylabel('$P_{\textrm{r}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
zlabel('$\omega$ (degrees)', 'fontsize', 20, 'Interpreter', 'LaTex')

% % plot P vs. s (colormap)
% subplot(1,2,1)
% surf(X,Y,s_mm')
% view(2)
% sbar = colorbar;
% ylabel(sbar, 's (mm)', 'fontsize', 20, 'Interpreter', 'LaTex')
% title('Linear Displacement', 'fontsize', 20)
% xlabel('$P_{\textit{l}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
% ylabel('$P_{\textrm{r}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
% zlabel('s (mm)', 'fontsize', 20, 'Interpreter', 'LaTex')
% 
% % plot P vs. w (colormap)
% subplot(1,2,2)
% surf(X,Y,w_deg')
% view(2)
% wbar = colorbar;
% ylabel(wbar, '$\omega$ (degrees)', 'fontsize', 20, 'Interpreter', 'LaTex')
% title('Rotational Displacement', 'fontsize', 20)
% xlabel('$P_{\textit{l}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
% ylabel('$P_{\textrm{r}}$ (kPa)', 'fontsize', 20, 'Interpreter', 'LaTex')
% zlabel('$\omega$ (degrees)', 'fontsize', 20, 'Interpreter', 'LaTex')
% 

