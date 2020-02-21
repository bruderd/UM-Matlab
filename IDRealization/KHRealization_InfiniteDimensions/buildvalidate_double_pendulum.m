%% simulate pendulum data to build model
h = 0.01;
y0 = [ 0; 0; 0; 0 ];
t0 = 0;
tfinal = 80;
input_t = t0:h:tfinal;
input_u =  -1 + 2.* rand( 2, length( input_t ) );

yout = ode1( @double_pendulum,t0,h,tfinal,y0,input_t,input_u,[] );

%% transfer function identification
Nout = size( yout, 2 );
Ninput = 2;
A = zeros( Nout, Ninput * Nout );

for i = 1:Nout
    A( i, ( Ninput * i ):-1:1 ) = reshape( input_u( :, 1: i ), 1, Ninput * i );
end

z = pinv( A ) * sin( yout( 1, : )' );

%% validation
input_u2 = -1 + 2.* rand( 2, length( input_t ) );

yout2 = ode1( @double_pendulum,t0,h,tfinal,y0,input_t,input_u2,[] );

%% input matrix
B = zeros( Nout, Ninput * Nout );

for i = 1:Nout
    B( i, ( Ninput * i ):-1:1 ) = reshape( input_u2( :, 1: i ), 1, Ninput * i );
end

figure; plot( input_t, abs( B * z - sin( yout2( 1, : ) )')./norm( sin( yout2( 1, : ) ) ) );
set(gca, 'YScale', 'log')
