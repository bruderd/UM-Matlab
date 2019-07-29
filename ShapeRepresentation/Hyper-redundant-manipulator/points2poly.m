function [coeffs, obs_matrix] = points2poly(degree, points, positions, orient)
%points2poly: Finds polynomial that best goes through a set of points.
% Polynomial starts at the origin, and its domain is [0,1].
% "Resting configuration" is along the yaxis (2d) or zaxis (3d)
%   degree - scalar, maximum degree of the polynomial
%   points - matrix, rows are xy(z) coordinates of points
%   positions - scaler [0,1], position of point on the arm.
%   orient - vector, orientation of the the end effector (complex number for 2d case, quaternion for 3d case)
%   coeffs - matrix, rows are the coefficients of the polynomial in each
%            coordinate where [a b c ...] -> ax^1 + bx^2 + cx^2 + ...
%   obs_matrix - matrix that converts from state vector to coefficients 

%% 2d case
if size(points,2) == 2
    if all( size(orient) ~= [1,2] )
        error('orientation for 2d system must be a complex number specified as [1x2] vector');
    end
    
    % generate virtual points to provide slope constraint at the base & end
    startpoint = [ 0 , 1e-2 ];
    endpoint = complex_mult( orient/norm(orient) , [ 0 , 1 ] )*1e-2 + points(end,:);
    points_supp = [0 , 0 ; startpoint ; points ; endpoint];
%     points_supp = points;   % remove the slope constraints
    
    % generate A matrix for least squares problem
    positions_supp = [ 0 , 1e-2 , positions , 1+1e-2 ];
%     positions_supp = positions;   % remove the slope constraints
    A = zeros( length(positions_supp) , degree );
    for i = 1 : degree
        A(:,i) = reshape( positions_supp , [] , 1) .^ i;
    end
    
    % separate x and y corrdinates of points
    points_x = points_supp(:,1);
    points_y = points_supp(:,2);
    
    % find polynomial coefficients
    obs_matrix = pinv(A);
    coeffs_vec_x = obs_matrix * points_x;
    coeffs_vec_y = obs_matrix * points_y;
    
    % make coeffs a matrix to where each row is coeffs for one dimension
    coeffs = [ coeffs_vec_x' ; coeffs_vec_y' ];
    
else
    error('points matrix must be nx2'); 
end

%% For generating plots (need to insert a breakpoint to use)
% 
% px = fliplr( [ 0 , coeffs(1,:) ] );
% py = fliplr( [ 0 , coeffs(2,:) ] );
% 
% for i = 1:101
%     pol_x(i) = polyval(px,0.01*(i-1));
%     pol_y(i) = polyval(py,0.01*(i-1));
%     s(i) = 0.01*(i-1);
% end

end