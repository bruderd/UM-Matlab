function axang = tform2axang( H )
%TFORM2AXANG Extract axis-angle rotation from homogeneous transformation
%   AXANG = TFORM2AXANG(H) retrieves the rotational component of the
%   3D homogeneous transformation, H, and returns it in an axis-angle
%   representation, AXANG. The translational components of H
%   will be ignored.
%   The input, H, is an 4-by-4-by-N matrix of N homogeneous transformations.
%   The output, AXANG, is an N-by-4 matrix of N axis-angle rotations.
%   The first three element of every row specify the rotation axis and
%   the last element defines the rotation angle (in radians).
%
%   Example:
%      % Convert a homogeneous transformation into an axis-angle vector
%      H = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1];
%      axang = tform2axang(H)
%
%   See also axang2tform

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateHomogeneousTransform(H, 'tform2axang', 'H');

% This is a two-step process.
% 1. Extract the rotation matrix from the homogeneous transform
R = tform2rotm(H);

% 2. Convert the rotation matrix into the axis-angle representation
axang = rotm2axang(R);

end

