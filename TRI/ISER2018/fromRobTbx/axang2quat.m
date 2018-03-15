function q = axang2quat(axang)
%AXANG2QUAT Convert axis-angle rotation representation to quaternion
%   Q = AXANG2QUAT(AXANG) converts a 3D rotation given in axis-angle form,
%   AXANG, to its unit quaternion representation, Q. AXANG is an N-by-4
%   matrix of N axis-angle rotations. The first three element of every
%   row specify the rotation axis and the last element defines the rotation
%   angle (in radians).
%   The output, Q, is an N-by-4 matrix containing N quaternions. Each quaternion 
%   is of the form q = [w x y z], with a scalar number as the first value.
%
%   Example:
%      % Convert a rotation from axis-angle to quaternion
%      axang = [0 1 0 pi/2];
%      q = axang2quat(axang)
%
%   See also quat2axang

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(axang, 'axang2quat', 'axang', ...
    'ncols', 4);

% For a single axis-angle vector [ax ay az t] the output quaternion
% q = [w x y z], can be computed as follows:
% w = cos(t/2)
% x = ax*sin(t/2)
% y = ay*sin(t/2)
% z = az*sin(t/2)

% Normalize the axis
v = robotics.internal.normalizeRows(axang(:,1:3));

% Create the quaternion
thetaHalf = axang(:,4)/2;
sinThetaHalf = sin(thetaHalf);
q = [cos(thetaHalf), v(:,1).*sinThetaHalf, v(:,2).*sinThetaHalf, v(:,3).*sinThetaHalf];

end
