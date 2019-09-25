function AB = matrix_mult(A,B)
%matrix_mult: Performs matrix multiplication using a for loop
%   Using this to see if Matlab's built in matrix multiplication yields a
%   different result.

if size(A,2) ~= size(B,1)
    error('Matrices to be multiplied are not consistent in dimension')
end

AB = zeros( size(A,1) , size(B,2) );
for i = 1 : size(A,1)
    for j = 1 : size(B,2)
        AB(i,j) = A(i,:) * B(:,j);
    end
end

end

