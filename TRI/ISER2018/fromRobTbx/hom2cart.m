function C = hom2cart(H)
%HOM2CART - Convert homogeneous coordinates to Cartesian coordinates
%   C = HOM2CART(H) converts a set of homogeneous points, H, into Cartesian
%   coordinates, C. H is of size N-by-K and contains N homogeneous points.
%   K needs to be greater or equal to 2.
%   The output, C, is an N-by-(K-1) matrix and contains N points in
%   Cartesian coordinates. Each row of C represents a point in
%   (K-1)-dimensional space.
%
%   Example:
%      % Convert two homogeneous points to 3D Cartesian points
%      h = [0.2785 0.9575 0.1576 0.5; 0.5469 0.9649 0.9706 0.5];
%      c = hom2cart(h)
%
%   See also CART2HOM.

%   Copyright 2014 The MathWorks, Inc.

%#codegen

robotics.internal.validation.validateNumericMatrix(H, 'hom2cart', 'H');

% Check minimum number of 2 columns in input
if size(H, 2) < 2
    error(message('robotics:validation:InvalidMinColumns', 'H'));
end

% Normalizes according to last column and removes it
d = size(H);
C = H(:,1:end-1) ./ repmat(H(:,end), 1, d(2)-1);

end

