function tr = tform2trvec( H )
%TFORM2TRVEC Extract translation vector from homogeneous transformation
%   T = TFORM2TRVEC(H) extracts the translation vector, T, from a 3D homogeneous
%   transformation, H, and returns it. The rotational components of H
%   will be ignored.
%   The input, H, is an 4-by-4-by-N matrix of N homogeneous transformations.
%   The output, T, is an N-by-3 matrix containing N translation vectors. Each
%   vector is of the form t = [x y z].
%
%   Example:
%      % Extract translation vector from homogeneous transformation
%      H = [1 0 0 0.5; 0 -1 0 5; 0 0 -1 -1.2; 0 0 0 1];
%      t = tform2trvec(H)
%
%   See also trvec2tform

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateHomogeneousTransform(H, 'tform2trvec', 'H');

% Also normalize by last element in matrix
t = H(1:end-1,end,:) ./ repmat(H(end,end,:), [3 1 1]);
tr = permute(t,[3 1 2]);

end

