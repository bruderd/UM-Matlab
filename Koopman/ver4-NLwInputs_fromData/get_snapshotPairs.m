function [ x, y ] = get_snapshotPairs(data, params)
%get_snapshotPairs: Generate a bunch of snapshot pairs from the data by
%   Detailed explanation goes here

% initialize output arrays
x = []; y = [];

for j = 1:length(data.x)-1
    x(j,:) = [ data.x(j,:), data.u(j,:) ];
    y(j,:) = [ data.x(j+1,:), data.u(j,:) ];
end

% % append snapshot pairs from this trial onto set of all pairs
% x = [x; xk];
% y = [y; yk];


end