h = 0.01;
y0 = [ 0; 0 ];
t0 = 0;
tfinal = 10;
input_t = t0:h:tfinal;
input_u = -2 + 4.* rand( length( input_t ), 1 );

yout = ode1( @pendulum,t0,h,tfinal,y0,input_t,input_u,[] );
