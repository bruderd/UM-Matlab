function [x] = matrix2vector(A)
% Converts matrix into a column vector by stacking columns on top of each
% other sequentially

    [m, n] = size(A);
    
    x = zeros(m*n , 1);
    
    for k = 1:n
        x((k-1)*m+1 : (k-1)*m+m , 1) = A(:,k);
    end
    
end