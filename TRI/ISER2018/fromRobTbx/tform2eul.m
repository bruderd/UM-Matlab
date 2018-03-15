function eul = tform2eul( H, varargin )
%TFORM2EUL Extract Euler angles from homogeneous transformation
%   EUL = TFORM2EUL(H) extracts the rotational component from a 3D homogeneous
%   transformation, H, and returns it as Euler angles, EUL. The translational
%   components of H will be ignored.
%   H is an 4-by-4-by-N array containing N homogeneous transformation
%   matrices. The output EUL is an N-by-3 array of Euler rotation angles. 
%   Rotation angles are in radians.
%
%   EUL = TFORM2EUL(H, SEQ) calculates the Euler angles EUL for homogeneous
%   transformation H, and a specified rotation sequence, SEQ.
%
%   The default rotation sequence is 'ZYX', where the order of rotation
%   angles is Z Axis Rotation, Y Axis Rotation, and X Axis Rotation.
%
%   The following rotation sequences, SEQ, are supported: 'ZYX' and 'ZYZ'.
%
%   Example:
%      % Calculate the Euler angles from a transformation matrix
%      % By default, the ZYX axis order will be used.
%      H = [1 0 0 0.5; 0 -1 0 5; 0 0 -1 -1.2; 0 0 0 1];
%      eul = tform2eul(H)
%
%      % Calculate H based on a ZYZ rotation
%      eulzyz = tform2eul(H, 'ZYZ')
%
%   See also eul2tform

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateHomogeneousTransform(H, 'tform2eul', 'H');
seq = robotics.internal.validation.validateEulerSequence(varargin{:});

% This is a two-step process.
% 1. Extract the rotation matrix from the homogeneous transform
R = tform2rotm(H);

% 2. Convert the rotation matrix to a set of euler angles
eul = rotm2eul(R, seq);

end

