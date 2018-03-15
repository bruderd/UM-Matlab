function H = quat2tform( q )
%QUAT2TFORM Convert quaternion to homogeneous transformation
%   H = QUAT2TFORM(Q) converts a unit quaternion, Q, into a homogeneous
%   transformation matrix, H. The input, Q, is an N-by-4 matrix containing N 
%   quaternions. Each quaternion represents a 3D rotation and is of the form 
%   q = [w x y z], with a scalar number as the first value. Each element 
%   of Q must be a real number.
%   The output, H, is an 4-by-4-by-N matrix of N homogeneous transformations.
%
%   Example:
%      % Convert a quaternion to homogeneous transform
%      q = [0.7071 0.7071 0 0];
%      H = quat2tform(q)
%
%   See also tform2quat

%   Copyright 2014 The MathWorks, Inc.

%#codegen

% This is a two-step process.
% 1. Convert the quaternion input into a rotation matrix
R = quat2rotm(q);

% 2. Convert the rotation matrix into a homogeneous transform
H = rotm2tform(R);

end

