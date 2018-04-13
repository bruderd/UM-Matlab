% Plot a single frame of markers on a 3D cartesian plot
% Tim Dorn
% 22nd Mar 2008
% 
% --------------------------------------------------------------------
% Usage: plotmarkers(data, faceCol*, outerCol*)
% --------------------------------------------------------------------
% 
% Inputs:   data = a vector of length 3n where n = number of markers
%           faceCol* = [R, G, B] of marker face colour
%           outerCol* = [R, G, B] of marker outline colour
% 
% --------------------------------------------------------------------

function plotmarkers(data, faceCol, outerCol)

usage = 'Usage: plotmarkers(data, faceCol*, outerCol*)';

% input error checking
if nargin == 1,
    faceCol = [0,0,1];
    outerCol = [1,0,0];
elseif nargin == 2,
    outerCol = [1,0,0];
elseif nargin ~= 3,
    disp(usage)
end

if mod(length(data), 3) ~= 0,
    fprintf(stderr, 'incorrect input vector dimensions');
    return
end

l = length(data);
x = []; y = []; z = [];

% get list of 3d points
for i = 1:l
    if mod(i,3) == 1, x = [x ; data(i)]; end
    if mod(i,3) == 2, y = [y ; data(i)]; end
    if mod(i,3) == 0, z = [z ; data(i)]; end
end

% set views
camUp = [-1, 0, 0 ; ...
         -0.6 0.95 0 ; ...
         0, 1, 0 ; ...
         0, 1, 0];

camPos = [681.561 17536.108 0 ; ...
         9369 2524 -3570 ; ...
         7000.885 981.347 0 ; ...
         730.726 985.614 17320.508];   
     

% plot markers & ground
for i = 1:4,
    subplot(2,2,i, 'CameraUpVector', camUp(i,:), ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'CameraPosition', camPos(i,:));
    hold on
    grid on
    
    % markers
    plot3(x, y, z, 'MarkerFaceColor', faceCol ,'Marker', 'o', ...
        'LineStyle', 'none', 'Color', outerCol);
 
    % ground (same XZ dimensions as my modified treadmill.vtp in OpenSim)
    xg = [-500 2000 2000 -500]';
    yg = [0 0 0 0]';
    zg = [-300 -300 300 300]';
    c = [1 1 1 1]';
    fill3(xg, yg, zg, c)
    
    size = 12;
    xlabel('X', 'FontSize', size)
    ylabel('Y', 'FontSize', size)
    zlabel('Z', 'FontSize', size) 
end

view3d rot
