function [x_scaled, y_scaled, scaleup] = scale_snapshotPairs( snapshotPairs , params )
%scale_snapshotPairs: Summary of this function goes here
%   Detailed explanation goes here

% maximum value for each of the states (doesn't include input), in a row vector
maxStates = max( abs( max(snapshotPairs.y(:,1:params.n), [],1) ) , abs( min(snapshotPairs.y(:,1:params.n), [],1) ));

% scale each state so that it never goes outside the range [-1, 1]
x_scaled = snapshotPairs.x ./ [maxStates , ones(1, params.p)];
y_scaled = snapshotPairs.y ./ [maxStates , ones(1, params.p)];

% create nxn diagonal matrix that will scale the states back up after the Koopman operator has been learned
scaleup = diag(maxStates);

end

