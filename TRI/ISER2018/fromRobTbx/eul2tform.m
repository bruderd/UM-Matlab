function H = eul2tform( eul, varargin )
%EUL2TFORM Convert Euler angles to homogeneous transformation
%   H = EUL2TFORM(EUL) converts a set of 3D Euler angles, EUL, into a homogeneous
%   transformation matrix, H. EUL is an N-by-3 matrix of Euler rotation
%   angles. The output H is an 4-by-4-by-N matrix of N homogeneous transformations.
%
%   H = EUL2TFORM(EUL, SEQ) converts 3D Euler angles to a homogeneous transformation.
%   The Euler angles are specified by the axis rotation sequence, SEQ.
%
%   The default rotation sequence is 'ZYX', where the order of rotation
%   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
%
%   The following rotation sequences, SEQ, are supported: 'ZYX' and 'ZYZ'.
%
%   Example:
%      % Calculate the transformation matrix for a set of Euler angles
%      % By default, the ZYX axis order will be used.
%      angles = [0 pi/2 0];
%      H = eul2tform(angles)
%
%      % Calculate H based on a ZYZ rotation
%      Hzyz = eul2tform(angles, 'ZYZ')
%
%   See also tform2eul

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(eul, 'eul2quat', 'eul', ...
    'ncols', 3);
seq = robotics.internal.validation.validateEulerSequence(varargin{:});

% This is a two-step process.
% 1. Convert the Euler angles into a rotation matrix
R = eul2rotm(eul, seq);

% 2. Convert the rotation matrix into a homogeneous transform
H = rotm2tform(R);

end


