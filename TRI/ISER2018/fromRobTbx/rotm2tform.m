function H = rotm2tform( R )
%ROTM2TFORM Convert rotation matrix to homogeneous transform
%   H = ROTM2TFORM(R) converts the 3D rotation matrix, R, into a homogeneous
%   transformation, H. H will have no translational components.
%   R is an 3-by-3-by-N matrix containing N rotation matrices.
%   Each rotation matrix has a size of 3-by-3 and is orthonormal.
%   The output, H, is an 4-by-4-by-N matrix of N homogeneous transformations.
%
%   Example:
%      % Convert a rotation matrix to a homogeneous transformation
%      R = [1 0 0 ; 0 -1 0; 0 0 -1]
%      H = rotm2tform(R)
%
%   See also tform2rotm

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

% Ortho-normality is not tested, since this validation is expensive
robotics.internal.validation.validateRotationMatrix(R, 'rotm2tform', 'R');

numMats = size(R,3);

% The rotational components of the homogeneous transformation matrix
% are located in elements H(1:3,1:3).
H = zeros(4,4,numMats,'like',R);
H(1:3,1:3,:) = R;
H(4,4,:) = ones(1,1,numMats,'like',R);

end

