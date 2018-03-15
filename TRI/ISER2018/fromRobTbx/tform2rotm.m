function R = tform2rotm( H )
%TFORM2ROTM Extract rotation matrix from homogeneous transformation
%   R = TFORM2ROTM(H) extracts the rotational component from a 3D homogeneous
%   transformation, H, and returns it as an orthonormal rotation matrix, R.
%   The translational components of H will be ignored.
%   The input, H, is an 4-by-4-by-N matrix of N homogeneous transformations.
%   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices.
%
%   Example:
%      % Convert a homogeneous transformation in a rotation matrix
%      H = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1];
%      R = tform2rotm(H)
%
%   See also rotm2tform

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateHomogeneousTransform(H, 'tform2rotm', 'H');

R = H(1:3,1:3,:);

end

