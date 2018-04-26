function Boxpoints = setBoxtestPoints( center, width, height )
%setMtestPoints: Creates a file of test points in the shape of simple
%square with a point in the middle
%   Detailed explanation goes here

% set the increments between points
dw = width / 2;
dh = height / 2;

%% Make first row zeros for calibration purposes

points(1,:) = zeros(1,2);

for i = 1:10    %repeat each point 10 times
   points(5*(i-1)+2,:) = center;
   points(5*(i-1)+3,:) = center + [dw, dh];  
   points(5*(i-1)+4,:) = center + [dw, -dh]; 
   points(5*(i-1)+5,:) = center + [-dw, -dh]; 
   points(5*(i-1)+6,:) = center + [-dw, dh];     
end

%% Set output (6 dimensional points)
Boxpoints = zeros(length(points(:,1)),6);
Boxpoints(:,3:4) = points;