function dy =  pendulum(t,y,input_t,input_u,param)

dy = zeros( 2, 1 );
dy( 1 ) = y( 2 );

idx =  find( abs( t - input_t ) < 1e-10 );

dy( 2 ) = input_u( idx ) - sin( y( 1 ) );