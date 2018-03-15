function H = axang2tform( axang )
%AXANG2TFORM Convert axis-angle rotation representation to homogeneous transform
%   H = AXANG2TFORM(AXANG) converts a 3D rotation given in axis-angle form,
%   AXANG, to a homogeneous transformation matrix, H. AXANG is an N-by-4
%   matrix of N axis-angle rotations. The first three elements of every
%   row specify the rotation axis and the last element defines the rotation
%   angle (in radians).
%   The output H is an 4-by-4-by-N matrix of N homogeneous transformations.
%
%   Example:
%      % Convert a rotation from axis-angle to homogeneous transform
%      axang = [0 1 0 pi/2];
%      H = axang2tform(axang)
%
%   See also tform2axang

%   Copyright 2014 The MathWorks, Inc.

%#codegen

% This is a two-step process.
% 1. Convert the axis-angle input into a rotation matrix
R = axang2rotm(axang);

% 2. Convert the rotation matrix into a homogeneous transform
H = rotm2tform(R);

end

