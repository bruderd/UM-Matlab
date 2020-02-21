function dy =  linear_system_discrete(t,y,input_t,input_u,param)

dy = zeros( 2, 1 );
idx =  find( abs( t - input_t ) < 1e-10 );

% dy( 1 ) = y( 2 );
% dy( 2 ) = input_u( idx ) - y( 1 ) + 0.1 * y( 2 );

A = [ 1 , 0.01 ; -0.01 , 1]; B = [0 ; 0.01];
dy = A * y + B * input_u(idx);