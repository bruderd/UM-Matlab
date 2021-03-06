function cost = costfun_Tonly( arg , A , B , data , sys )
%costfun: Cost function for nonlinear optimization problem
%   arg - optimization decision variable
%   sys - sysid object
%   data - struct containing lifted snapshot pairs and inputs (Px, Py, u)
%   params - system parameters such as sizes of things (N, m, k, etc)

params = sys.params;
model = sys.model;
N = params.N;
m = params.m;

% convert decision variable into matrices
T = reshape( arg( 1 : N^2 ) , [N,N] );

% error for each snapshot pair
err = ( data.Px * T' * A' + data.u * B' - data.Py * T' );
L2err = vecnorm( err , 2 , 2 );
L1err = sum( abs( err ) , 2 );

% cost is the sum of l2 error over all points
cost = sum( L2err );

end

