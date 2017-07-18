% /simpleParallel/vfplot_forcebalance.m
%   Plots the vector fields for an axis aligned parallel configuration of two 1-fiber FREEs. Uses the force balance
%   equations for force/torque rather than the volume/energy derivations
%
%   Must specify parameters in setParms_2parallel.m before running this script.

[Gama, R, L, B, N] = deal(params.Gama_l, params.R_l, params.L_l, params.B_l, params.N_l);
[Gama_r, R_r, L_r, B_r, N_r] = deal(params.Gama_r, params.R_r, params.L_r, params.B_r, params.N_r);
[Kf, Km] = deal(params.kelast_l(1), params.kelast_l(2));
[F_load, M_load] = deal(params.load(1), params.load(2));
Ptest_l = params.Ptest_l;
Ptest_r = params.Ptest_r;

% scaling factors to make vector plots more legible
scaleF = 10e-4;
scaleM = 2*10e1;

%% Make vector field for pressure induced forces

smin = -(B-L)+(L/10);
smax = (B-L)-(L/10);
wmin = -2*pi;
wmax = 2*pi;
res_s = (smax-smin)/8;   % resolution in s
res_w = (wmax-wmin)/8;    % resolution in w

[s,w] = meshgrid(smin:res_s:smax, wmin:res_w:wmax);

% force and torque in terms of P, s, w

% forces of left FREE
dVds_l = pi*(B^2 - 3.*(L+s).^2) ./ (2*pi*N + w).^2;
dVdw_l = 2*pi*((L+s) .* ((L+s).^2 - B^2)) ./ (2*pi*N + w).^3;
Force_l = Ptest_l * dVds_l;
Moment_l = Ptest_l * dVdw_l;

% forces of right FREE
dVds_r = -pi*(B^2 - 3.*(L-s).^2) ./ (2*pi*N - w).^2;           % right FREE
dVdw_r = -2*pi*((L-s) .* ((L-s).^2 - B^2)) ./ (2*pi*N - w).^3; % right FREE
Force_r = Ptest_r * dVds_r;
Moment_r = Ptest_r * dVdw_r;

% scale for better plotting
Force_l_plot = Force_l * scaleF;
Moment_l_plot = Moment_l * scaleM;

Force_r_plot = Force_r * scaleF;
Moment_r_plot = Moment_r * scaleM;


%% Make vector field for f_elast
Felast = Kf * s;
Melast = Km * w;

% scale for better plotting
Felast_plot = Felast * scaleF;
Melast_plot = Melast * scaleM;


%% Make vector field for f_load
Fload = F_load * ones(size(s));
Mload = M_load * ones(size(w));

% scale for better plotting
Fload_plot = Fload * scaleF;
Mload_plot = Mload * scaleM;


%% Make vector field for net force
netF = Force_l + Force_r + Felast + Fload;
netM = Moment_l + Moment_r + Melast + Mload;

% scale for better plotting
netF_plot = netF * scaleF;
netM_plot = netM * scaleM;


%% Plot everything in one figure

figure
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)

subplot(2,2,1)
forcefield = quiver(s,w,Force_l_plot,Moment_l_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{l}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$ (m)','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$ (rad)','Interpreter','Latex', 'FontSize',25)
forcefield.MaxHeadSize = 0.005;
forcefield.Marker = '.';
forcefield.ShowArrowHead = 'off';
forcefield.Color = [160 180 180] ./ 255;

starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,Force_l_plot,Moment_l_plot,startx,starty)

subplot(2,2,2)
forcefield = quiver(s,w,Force_r_plot,Moment_r_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{r}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$ (m)','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$ (rad)','Interpreter','Latex', 'FontSize',25)
forcefield.MaxHeadSize = 0.005;
forcefield.Marker = '.';
forcefield.ShowArrowHead = 'off';
forcefield.Color = [160 180 180] ./ 255;

starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,Force_r_plot,Moment_r_plot,startx,starty)

subplot(2,2,3)
elastfield = quiver(s,w,Felast_plot,Melast_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{elastomer}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$ (m)','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$ (rad)','Interpreter','Latex', 'FontSize',25)
elastfield.MaxHeadSize = 0.005;
elastfield.Marker = '.';
elastfield.ShowArrowHead = 'off';
elastfield.Color = [160 180 180] ./ 255;

% subplot(2,2,3)
% netfield = quiver(s,w,Fload_plot,Mload_plot);
% axis([smin smax wmin wmax])
% title('$\mathbf{f}_{\textrm{load}}$','Interpreter','Latex', 'FontSize',25)
% xlabel('$s$ (m)','Interpreter','Latex', 'FontSize',25)
% ylabel('$\omega$ (rad)','Interpreter','Latex', 'FontSize',25)
% netfield.MaxHeadSize = 0.005;
% netfield.Marker = '.';
% netfield.ShowArrowHead = 'off';
% netfield.Color = [160 180 180] ./ 255;

subplot(2,2,4)
netfield = quiver(s,w,netF_plot,netM_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{net}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$ (m)','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$ (rad)','Interpreter','Latex', 'FontSize',25)
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

% plot trajectories from s = 0
starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,netF_plot,netM_plot,startx,starty)

% plot trajectories from w = 0
startx = smin:res_s:smax;
starty = zeros(size(startx));
streamline(s,w,netF_plot,netM_plot,startx,starty)

% % plot trajectories from diagonal
% startx = smin:res_s:smax;
% starty = wmin:res_w:wmax;
% streamline(s,w,netF_plot,netM_plot,startx,starty)
% 
% % plot trajectories from other diagonal
% startx = smin:res_s:smax;
% starty = -(wmin:res_w:wmax);
% streamline(s,w,netF_plot,netM_plot,startx,starty)