function testPoints = setSometestPoints( z, phi, params )
%setSometestPoint: Creates a file of test points based on input
%   z and phi must be row vectors of same size. They specify points


len = length(z);

%% Make first row zeros for calibration purposes

points(1,:) = zeros(1,2);

for i = 1:10    %repeat each point 10 times
   points(len*(i-1)+1 : len*i,:) = [z' , phi'];   
end

%% Set output (6 dimensional points)
testPoints = zeros(length(points(:,1))+1,6);    % ensure first row is zeros
testPoints(2:end,3:4) = points;

% testPoints(1:end, 5) = params.x5_offset;    % add offset
% testPoints(1:end, 6) = params.x6_offset;    % add offset