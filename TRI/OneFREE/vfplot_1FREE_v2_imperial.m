% /SimpleParallel/vfplot_forcebalance.m
%   Plots the vector fields for a 1-fiber FREE. Uses the force balance
%   equations for force/torque rather than the volume/energy derivations

[Gama, R, L, B, N] = deal(params.Gama, params.R, params.L, params.B, params.N);
[Kf, Km] = deal(params.kelast(1), params.kelast(2));
[F_load, M_load] = deal(params.load(1), params.load(2));
Ptest = params.Ptest;

% scaling factors to make plots more legible
scaleF = 10e0;
scaleM = 5*10e1;

%% Make vector field for pressure induced forces

smin = -(B-L)+(L/10);
smax = (B-L)-(L/10);
wmin = -2*pi;
wmax = 2*pi;
res_s = (smax-smin)/8;   % resolution in s
res_w = (wmax-wmin)/8;    % resolution in w

[s,w] = meshgrid(smin:res_s:smax, wmin:res_w:wmax);

% force and torque in terms of P, s, w
dVds = pi*(B^2 - 3.*(L+s).^2) ./ (2*pi*N + w).^2;
dVdw = 2*pi*((L+s) .* ((L+s).^2 - B^2)) ./ (2*pi*N + w).^3;
Force = Ptest * dVds;
Moment = Ptest * dVdw;

% scale for better plotting
Force_plot = Force * scaleF;
Moment_plot = Moment  * scaleM;


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


%% Make net force vector field
netF = Force + Felast + Fload;
netM = Moment + Melast + Mload;

% scale for better plotting
netF_plot = netF * scaleF;
netM_plot = netM * scaleM;


%% Plot everything in one figure

figure
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)

subplot(2,2,1)
forcefield = quiver(s,w,Force_plot,Moment_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{fiber}} + \mathbf{f}_{\textrm{endcap}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
forcefield.MaxHeadSize = 0.005;
forcefield.Marker = '.';
forcefield.ShowArrowHead = 'off';
forcefield.Color = [160 180 180] ./ 255;

starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,Force_plot,Moment_plot,startx,starty)

subplot(2,2,2)
elastfield = quiver(s,w,Felast_plot,Melast_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{elastomer}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
elastfield.MaxHeadSize = 0.005;
elastfield.Marker = '.';
elastfield.ShowArrowHead = 'off';
elastfield.Color = [160 180 180] ./ 255;

subplot(2,2,3)
netfield = quiver(s,w,Fload_plot,Mload_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{load}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

subplot(2,2,4)
netfield = quiver(s,w,netF_plot,netM_plot);
axis([smin smax wmin wmax])
title('$\mathbf{f} = \mathbf{f}_{\textrm{fiber}} + \mathbf{f}_{\textrm{endcap}} + \mathbf{f}_{\textrm{elastomer}} + \mathbf{f}_{\textrm{load}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

% plot trajectories from diagonal
startx = smin:res_s:smax;
starty = wmin:res_w:wmax;
streamline(s,w,netF_plot,netM_plot,startx,starty)

% plot trajectories from other diagonal
startx = smin:res_s:smax;
starty = -(wmin:res_w:wmax);
streamline(s,w,netF_plot,netM_plot,startx,starty)

% starty = wmin:res_w:wmax;
% startx = zeros(size(starty));
% streamline(s,w,netF,netM,startx,starty)
% 
% startx = smin:res_s:smax;
% starty = ones(size(startx));
% streamline(s,w,netF,netM,startx,starty)