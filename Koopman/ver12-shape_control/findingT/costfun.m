function cost = costfun( arg , data , params )
%costfun: Cost function for nonlinear optimization problem
%   arg - optimization decision variable
%   data - struct containing lifted snapshot pairs and inputs (Px, Py, u)
%   params - system parameters such as sizes of things (N, m, k, etc)

N = params.N;
m = params.m;

% convert decision variable into matrices
T = reshape( arg( 1 : N^2 ) , [N,N] );
A = reshape( arg( N^2 + 1 : 2*N^2 ) , [N,N] );
B = reshape( arg( 2*N^2 + 1 : end ) , [N,m] );

% error for each snapshot pair
err = data.Px * T' * A' + data.u * B' - data.Py * T';
L2err = vecnorm( err , 2 , 2 );
L1err = sum( abs( err ) , 2 );

% cost is the sum of l2 (or l1) error over all points
cost = sum( L2err );

end

