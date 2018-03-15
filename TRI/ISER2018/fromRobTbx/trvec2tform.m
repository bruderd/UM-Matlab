function H = trvec2tform( t )
%TRVEC2TFORM Convert translation vector to homogeneous transformation
%   H = TRVEC2TFORM(T) converts the Cartesian representation of a 3D translation 
%   vector, T, into the corresponding homogeneous transformation, H.
%   The input, T, is an N-by-3 matrix containing N translation vectors. Each
%   vector is of the form t = [x y z]. The output, H, is an 4-by-4-by-N matrix 
%   of N homogeneous transformations.
%
%   Example:
%      % Create homogeneous transformation from translation vector
%      t = [0.5 6 100];
%      H = trvec2tform(t)
%
%   See also tform2trvec

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(t, 'trvec2tform', 't', ...
    'ncols', 3);

numTransl = size(t, 1);

H = repmat(eye(4,'like',t),[1,1,numTransl]);
H(1:end-1,end,:) = t.';

end

