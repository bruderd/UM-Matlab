function dy =  double_pendulum(t,y,input_t,input_u,param)

dy = zeros( 4, 1 );
dy( 1 ) = y( 2 );
dy( 3 ) = y( 4 );
idx =  find( abs( t - input_t ) < 1e-10 );

m1 = 1;
l1 = 1/2;
I1 = 1;
m2 = 1;
l2 = 1/2;
I2 = 1;
L = 1;
g = 1;

c1 = m1*l1^2/2 + I1/2 + m2*L^2/2;
c2 = m2*l2^2/2 + I2/2;
c3 = m2*L*l2;
c4 = g*(m1*l1 + m2*L);
c5 = g*m2*l2;

dy( 2 ) = ( 2 * c2 * c4 * sin( y( 1 ) ) + c3^2 * y( 2 )^2 * sin( y( 1 ) - y( 3 ) ) * cos( y( 1 ) - y( 3 ) ) + ...
    2 * c2 * c3 * y( 4 )^2 * sin( y( 1 ) - y( 3 ) ) - c3 * c5 * cos( y( 1 ) - y( 3 ) ) * sin( y( 3 ) ) )/ ...
    ( c3^2 * ( cos( y( 1 ) - y( 3 ) ) )^2 - 4 * c1 * c2 ) + input_u( 1, idx );


dy( 4 ) = ( 2 * c1 * c5 * sin( y( 3 ) ) - c3^2 * y( 4 )^2 * sin( y( 1 )- y( 3 ) ) * cos( y( 1 ) - y( 3 ) ) - ...
    2 * c1 * c3 * y( 2 )^2 * sin( y( 1 ) - y( 3 ) ) - c3 * c4 * cos( y( 1 ) - y( 3 ) ) * sin( y( 1 ) ) )/ ...
    ( c3^2 * ( cos( y( 1 ) - y( 3 ) ) )^2 - 4 * c1 * c2 ) + input_u( 2, idx );