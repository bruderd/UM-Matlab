function vf_visualization( csv_linear, csv_constitutive, res, xrange, yrange )
%vf_visualization
%   Reads in csv files (vf1, vf2, ...) that describe vector fields then
%   plots them both in the same figure

if nargin == 2
    res = 11;   % number of vector points in the final plot
    xrange = [-0.25,0.25];    % [xmin, xmax]
    yrange = [-1,1]; % [ymin, ymax]
elseif nargin == 3
    xrange = [-0.25,0.25];    % [xmin, xmax]
    yrange = [-1,1]; % [ymin, ymax]    
end

% Read in the relaxed FREE parameters (Gama, L, R) and the pressure (P)
params1 = csvread(csv_linear, 0, 0, 'A1..D1');
params2 = csvread(csv_constitutive, 0, 0, 'A1..D1');

% Make sure the both vector fields are for the same FREE 
if params1 ~= params2
    disp('Error: You are attempting to compare vector fields for FREEs with different parameters.');
    return;
end

% extract FREE parameters
[Gama, L, R, P] = deal(params1(1), params1(2), params1(3), params1(4));

% Read in the vector fields
vf1 = csvread(csv_linear, 1, 0);
vf2 = csvread(csv_constitutive, 1, 0);
[x1, y1, u1, v1] = deal( vf1(:,1), vf1(:,2), vf1(:,3), vf1(:,4) );
[x2, y2, u2, v2] = deal( vf2(:,1), vf2(:,2), vf2(:,3), vf2(:,4) );

% Change units of s and w such that s: % of L, w: rotations
s1 = x1 * (1/L);
s2 = x2 * (1/L);
w1 = y1 * (1/(2*pi));
w2 = y2 * (1/(2*pi));

% Scale F and M for better plotting
fscale = 3e-4; % changes apparent width of arrows
mscale = 5e-1; % changes apparent height of arrows
f1 = u1 * fscale;
f2 = u2 * fscale;
m1 = v1 * mscale;
m2 = v2 * mscale;


% Define the points on the vector field to be plotted
X = linspace(xrange(1), xrange(2), res);
Y = linspace(yrange(1), yrange(2), res);
[S,W] = meshgrid(X,Y);

% interpolate to find force vectors at the points to be plotted
F1 = griddata(s1, w1, f1, S,W);
M1 = griddata(s1, w1, m1, S,W);

F2 = griddata(s2, w2, f2, S,W);
M2 = griddata(s2, w2, m2, S,W);


% create figure
vfplot = figure;
set(vfplot, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
title(['$\Gamma =$ ' num2str(Gama) ' (deg), $L =$ ' num2str(L*1e2) ' (cm), $R = $' num2str(R*1e2) ' (cm)'],'Interpreter','Latex', 'FontSize',25)
xlabel({'$s$ ($\%$ of L)' ; '$F$ ($3 \times 10^{-4}$ N)'},'Interpreter','Latex', 'FontSize',25)
ylabel({'$w$ (rev)' ; '$M$ ($5 \times 10^{-1}$ Nm)'},'Interpreter','Latex', 'FontSize',25)
axis([xrange(1) xrange(2) yrange(1) yrange(2)]);
line(xrange, [0,0]);    % draw x-axis
line([0,0], yrange);    % draw y-axis
grid on;
hold on
quiv1 = quiver(S,W, F1, M1, 'Color', 'red', 'Linewidth', 1, 'MaxHeadSize', 5e-2, 'AutoScale', 'off', 'ShowArrowHead', 'off', 'Marker', '.');
quiv2 = quiver(S,W, F2, M2, 'Color', 'blue', 'Linewidth', 1, 'MaxHeadSize', 5e-2, 'AutoScale', 'off', 'ShowArrowHead', 'off', 'Marker', '.');
hold off
% legend([quiv1, quiv2], {'Linear Model', 'Constitutive Model'}, 'Fontsize', 24, 'Location', 'eastoutside');


%% Add arrowheads that look good (need the arrow3 function to work)

hold on
U1 = reshape(quiv1.UData, [numel(quiv1.UData),1]);
V1 = reshape(quiv1.VData, [numel(quiv1.VData),1]);
U2 = reshape(quiv2.UData, [numel(quiv2.UData),1]);
V2 = reshape(quiv2.VData, [numel(quiv2.VData),1]);
X = reshape(quiv1.XData, [numel(quiv1.XData),1]);
Y = reshape(quiv1.YData, [numel(quiv1.YData),1]);

arrow3( [X,Y], [X+U1,Y+V1] , 'r', 0.5, 0.5, [], 1, []);
arrow3( [X,Y], [X+U2,Y+V2] , 'b', 0.5, 0.5, [], 1, []);
hold off

%% Draw legend for plot
legend([quiv1, quiv2], {'Linear Model', 'Constitutive Model'}, 'Fontsize', 24, 'Location', 'eastoutside');

end
