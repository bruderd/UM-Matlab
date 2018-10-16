function W = calc_W( L , dataPoints , params)
%calc_W: Calculates the coefficient matrix W that satisfies xdot = W*psi(x,u)
%   Detailed explanation goes here

n = params.n;       % dimension of state, x
N = params.N;       % length of the basis
K = params.K; % size(dataPoints,1);     % total number of datapoints

Ldiag = kron( ones(K,1) , L');    % diagonally stack the transpose of L

% evaluate the basis jacobian at each point and stack the result
dpsi_dx = zeros(K*N, n);
for i = 1 : K
    x = dataPoints( i , 1:n )';
    u = dataPoints( i , (n+1):end )';
    dpsi_dx( (i-1)*N+1 : i*N , : ) =  jacobianLift(x,u);
end

W = dpsi_dx \ Ldiag;

end

