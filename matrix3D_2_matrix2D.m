function [x] = matrix3D_2_matrix2D(A)
% Converts 3D matrix into a 2D matrix by stacking rows on top of each
% other sequentially, and turning the 3rd dimension into the column
% dimension.

    [m, n, p] = size(A);
    
    x = zeros(m*n , p);
    
    for k = 1:m
        x((k-1)*n+1 : (k-1)*n+n , :) = A(k,:,:);
    end

    x = x';
end