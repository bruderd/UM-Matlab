% /SimpleParallel/vfplot.m
%   Plots the vector fields for a 1-fiber FREE

[Gama, R, L, B, N] = deal(params.Gama, params.R, params.L, params.B, params.N);
[Kf, Km] = deal(params.kelast(1), params.kelast(2));
[F_load, M_load] = deal(params.load(1), params.load(2));
Ptest = params.Ptest;

%% Define Jacobain

smin = -L;
smax = L;
wmin = -2*pi;
wmax = 2*pi;
% res_s = 0.002;   % resolution in s
% res_w = 0.4;    % resolution in w
res_s = L/4;   % resolution in s
res_w = pi/4;    % resolution in w

[s,w] = meshgrid(smin:res_s:smax, wmin:res_w:wmax);
dVds = pi*(B^2 - 3.*(L+s).^2) ./ (2*pi*N + w).^2;
dVdw = 2*pi*((L+s) .* ((L+s).^2 - B^2)) ./ (2*pi*N + w).^3;

pdVds = (8e-1)*dVds;
pdVdw = (8e2)*dVdw;

% multiply by pressure to get force and torque
Force = Ptest * dVds;
Torque = Ptest *dVdw;

% plot vector field
figure
forcefield = quiver(s,w,Force,Torque);
axis([smin smax wmin wmax])
forcefield.MaxHeadSize = 0.005;
forcefield.Marker = '.';
forcefield.ShowArrowHead = 'off';
forcefield.Color = [160 180 180] ./ 255;

% plot trajectories from w = 0
starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,Force,Torque,startx,starty)

% plot trajectories from s = 0
startx = smin:res_s:smax;
starty = 3*ones(size(startx));
streamline(s,w,Force,Torque,startx,starty)

% % find where the Jacobian is "zero"
% locked_s = abs(dVds) < 1e-5;
% locked_w = abs(dVdw) < 1e-5;
% locked = locked_s & locked_w;
% 
% figure
% plotmatrix(s(locked), w(locked))
% axis([smin smax wmin wmax])

%% Make vector field for f_elast
Felast = Kf * s;
Melast = Km * w;

figure
elastfield = quiver(s,w,Felast,Melast);
axis([smin smax wmin wmax])
elastfield.MaxHeadSize = 0.005;
elastfield.Marker = '.';
elastfield.ShowArrowHead = 'off';

%% Make vector field for f_load
Fload = F_load * ones(size(s));
Mload = M_load * ones(size(w));

figure
elastfield = quiver(s,w,Fload,Mload);
axis([smin smax wmin wmax])
elastfield.MaxHeadSize = 0.005;
elastfield.Marker = '.';
elastfield.ShowArrowHead = 'off';


%% Make net force vector field
netF = Ptest*dVds + Felast + Fload;
netM = Ptest*dVdw + Melast + Mload;

figure
netfield = quiver(s,w,netF,netM);
axis([smin smax wmin wmax])
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,netF,netM,startx,starty)


%find where the net force is "zero"
netF0 = abs(netF) < 1e0;
netM0 = abs(netM) < 1e-2;
net0 = netF0 & netM0;

figure
plotmatrix(s(net0), w(net0))
axis([smin smax wmin wmax])

%% Plot everything in one figure

figure
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',30)

subplot(2,2,1)
forcefield = quiver(s,w,Force,Torque);
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
streamline(s,w,Force,Torque,startx,starty)

subplot(2,2,2)
elastfield = quiver(s,w,Felast,Melast);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{elastomer}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
elastfield.MaxHeadSize = 0.005;
elastfield.Marker = '.';
elastfield.ShowArrowHead = 'off';
elastfield.Color = [160 180 180] ./ 255;

subplot(2,2,3)
netfield = quiver(s,w,Fload,Mload);
axis([smin smax wmin wmax])
title('$\mathbf{f}_{\textrm{load}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

subplot(2,2,4)
netfield = quiver(s,w,netF,netM);
axis([smin smax wmin wmax])
title('$\mathbf{f} = \mathbf{f}_{\textrm{fiber}} + \mathbf{f}_{\textrm{endcap}} + \mathbf{f}_{\textrm{elastomer}} + \mathbf{f}_{\textrm{load}}$','Interpreter','Latex', 'FontSize',25)
xlabel('$s$','Interpreter','Latex', 'FontSize',25)
ylabel('$\omega$','Interpreter','Latex', 'FontSize',25)
netfield.MaxHeadSize = 0.005;
netfield.Marker = '.';
netfield.ShowArrowHead = 'off';
netfield.Color = [160 180 180] ./ 255;

starty = wmin:res_w:wmax;
startx = zeros(size(starty));
streamline(s,w,netF,netM,startx,starty)

startx = smin:res_s:smax;
starty = ones(size(startx));
streamline(s,w,netF,netM,startx,starty)
















