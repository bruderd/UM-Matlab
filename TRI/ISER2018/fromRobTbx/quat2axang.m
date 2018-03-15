function axang = quat2axang(q)
%QUAT2AXANG Convert quaternion to axis-angle rotation representation
%   AXANG = QUAT2AXANG(Q) converts a unit quaternion, Q, into the equivalent
%   axis-angle representation, AXANG. The input, Q, is an N-by-4 matrix
%   containing N quaternions. Each quaternion represents a 3D rotation
%   and is of the form q = [w x y z], with a scalar number as the first
%   value. Each element of Q must be a real number.
%   AXANG is an N-by-4 matrix of N axis-angle rotations. The first
%   three element of every row specify the rotation axis and the last
%   element defines the angle (in radians, in the interval [-pi pi]).
%
%   Example:
%      % Convert a quaternion to axis-angle
%      q = [0.7071 0.7071 0 0];
%      axang = quat2axang(q)
%
%   See also axang2quat

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(q, 'quat2axang', 'q', ...
    'ncols', 4);

% Normalize the quaternions
q = robotics.internal.normalizeRows(q);

% Normalize and generate the rotation vector and angle sequence
% For a single quaternion q = [w x y z], the formulas are as follows:
% (axis) v = [x y z] / norm([x y z]);
% (angle) theta = 2 * acos(w)
v = robotics.internal.normalizeRows(q(:,2:4));
theta = robotics.internal.wrapToPi(2*acos(q(:,1)));

% Handle the zero rotation degenerate case
% Can return an arbitrary axis, but fix on z-axis rotation
zeroIdx = abs(theta) < 10 * eps(class(q));

numZeroIdx = sum(zeroIdx);
assert(numZeroIdx <= length(zeroIdx));

if any(zeroIdx)
    v(zeroIdx,:) = repmat([0 0 1], numZeroIdx, 1);
    theta(zeroIdx) = 0;
end

axang = cat(2, v, theta);

end

