function [A,N,B,C] = get_realization_from_H( H , m , n )
%get_realization_from_H: Extracts a bilinear realization from the Hankel
%matrix of the form: x[k+1] = Ax[k] + Nx[k]u[k] + Bu[k] , y[k] = Cx[k]
%   
% Inputs:
%   H - Hankel matrix
%   m - "block" height of hankel matrix
%   n - "block" width of hankel matrix
%
% Note: SISO systems only.

% Remove last "block" dimension from width and height
Hs = H( 1 : 2^(m-1)-1 , 1 : 2^(n-1)-1 ); 

% Take SVD of Hankel matrix
[U,S,V] = svd( Hs , 'econ' );
largesvs = find( diag(S) > 1e-4 );
Q = U(:,largesvs) * S(largesvs,largesvs);    % notation to match paper
P = V(:,largesvs)';     % notation to match paper
% Q = U*S;    % notation to match paper
% P = V';     % notation to match paper

% Partition the Hankel matrix (see paper) (H1, H2, and Hs should have same dimensions)
H1 = zeros( size(Hs) ); % preallocate
H2 = zeros( size(Hs) ); % preallocate
for i = 1 : n - 1
    
    if i == 1
        col_index = 1;
        one_index = 2;
        two_index = 3;
    else
        col_index = 2^(i-1) : 2^i - 1;
        one_index = 2^i : 2^i + 2^(i-1) - 1;
        two_index = 2^i + 2^(i-1) : 2^(i+1) - 1;
    end
    
    H1( : , col_index ) = H( 1 : 2^(m-1)-1 , one_index );
    H2( : , col_index ) = H( 1 : 2^(m-1)-1 , two_index );
end

% Assign outputs
A = pinv( Q' * Q ) * Q' * H1 * P' * pinv( P * P' );
N = pinv( Q' * Q ) * Q' * H2 * P' * pinv( P * P' );
B = P(:,1);
C = Q(1,:);

end