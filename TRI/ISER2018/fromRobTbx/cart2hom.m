function H = cart2hom(C)
%CART2HOM - Convert Cartesian coordinates to homogeneous coordinates
%   H = CART2HOM(C) converts a set of points in Cartesian coordinates, C,
%   into homogeneous points, H. C is an N-by-K matrix and contains N
%   points in Cartesian coordinates. Each row of C represents one point
%   in K-dimensional space. K needs to be greater or equal to 1.
%   The output, H, is an N-by-(K+1) matrix containing N homogeneous
%   points.
%
%   Example:
%      % Convert two 3D Cartesian points to homogeneous coordinates
%      c = [0.8147 0.1270 0.6324; 0.9058 0.9134 0.0975];
%      h = cart2hom(c)
%
%   See also HOM2CART

%   Copyright 2014-2015 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(C, 'cart2hom', 'C');

nrows = size(C,1);
unitvec = ones(nrows,1, 'like', C);
H = cat(2,C, unitvec);

end

