function [x_scaled, y_scaled, L_scale, R_scale] = scale_snapshotPairs( snapshotPairs , params )
%scale_snapshotPairs: Summary of this function goes here
%   Detailed explanation goes here

% maximum value for each of the states (doesn't include input), in a row vector
maxStates = max( abs( max(snapshotPairs.y(:,1:params.n), [],1) ) , abs( min(snapshotPairs.y(:,1:params.n), [],1) ));

% scale each state so that it never goes outside the range [-1, 1]
x_scaled = snapshotPairs.x ./ [maxStates , 10 * ones(1, params.p)];
y_scaled = snapshotPairs.y ./ [maxStates , 10 * ones(1, params.p)];

% create nxn diagonal matrix that will scale left side of the coefficientmatrix
L_scale = diag(maxStates);

% create NxN diagonal matrix that scales the right side of coefficent matrix
maxStates_inv = 1 ./ maxStates;
scaledBasis = polyLift(maxStates_inv' , 0.1 * ones(params.p, 1));
R_scale = diag(scaledBasis);

end

