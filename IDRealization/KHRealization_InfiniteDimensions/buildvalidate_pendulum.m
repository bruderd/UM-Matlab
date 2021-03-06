%% simulate pendulum data to build model
h = 0.01;
slen = 100;
y0 = [ 0 ; 0 ];
t0 = 0;
tfinal = 300;
input_t = t0:h:tfinal;
utrain = 1 .* rand( ( length( input_t ) - 1 ) / slen , 1 );
input_u = zeros( slen * ( length(utrain)+1 ) , 1 );
for i = 2 : length(utrain)
    input_u( slen*(i-1)+1 : slen*i , : ) = utrain( i );
end
% input_u = [ 0 ; interleave2( utrain , utrain , utrain , utrain , utrain , utrain , utrain , utrain , 'row' ) ];

yout = ode1( @pendulum,t0,h,tfinal,y0,input_t,input_u,[] );
% yout = ode1_discrete( @linear_system_discrete,t0,h,tfinal,y0,input_t,input_u,[] );

% make output weirder
yout( 1, : ) = 1 - cos( yout( 1, : ) );
yout( 2, : ) = sin( yout( 1, : ) ) ;

%% transfer function identification
Nout = size( yout, 2 );
A = zeros( Nout );

for i = 1:Nout
    A( i, i:-1:1 ) = input_u(1:i);
end

z = pinv( A ) * yout( 1, : )';  % only finds tf for first output y1

%% validation
y02 = [ 0 ; 0 ];
tfinal2 = 10;
input_t2 = t0:h:tfinal2;
uval = 0.5 .* rand( ( length( input_t2 ) - 1 ) / slen , 1 );
input_u2 = zeros( slen * ( length(uval)+1 ) , 1 );
for i = 2 : length(uval)
    input_u2( slen*(i-1)+1 : slen*i , : ) = uval( i );
end
% input_u2 = [ 0 ; interleave2( uval , uval , uval , uval , uval , uval , uval , uval , 'row' ) ];
% input_u2 = 0.5 .* ones( length( input_t2 ) , 1 );
Nout2 = length( input_t2 );

yout2 = ode1( @pendulum,t0,h,tfinal2,y02,input_t2,input_u2,[] );
% yout2 = ode1_discrete( @linear_system_discrete,t0,h,tfinal2,y02,input_t2,input_u2,[] );

% make output weirder
yout2( 1, : ) = 1 - cos( yout2( 1, : ) );
yout2( 2, : ) = sin( yout2( 1, : ) ) ;

%% input matrix
B = zeros( Nout2 );

for i = 1:Nout2
    B( i, i:-1:1 ) = input_u2( 1: i );
end

% plot real response verses predicted response
figure;
hold on;
plot( input_t2 , yout2( 1 , : )' , 'LineWidth' , 2);
plot( input_t2 , B * z(1:Nout2) , 'LineWidth' , 2);
hold off;

% plot error
figure; plot( input_t2, abs( B * z(1:Nout2) - yout2( 1 , : )')./norm( yout2(1 , : ) ) );
set(gca, 'YScale', 'log')

%% CONSTRUCT HO-KALMAN REALIZATION (need to fix this)

% build block Hankel matrix from impulse response
H = hankel( z(1:floor(Nout2/2)) , z(floor(Nout2/2):Nout2) );
% Hrank = rank(H);
Hrank = 2;

% compute svd of the Hankel matrix
[ U , S , V ] = svd( H );

% estimate the controllability and observability matrices
% Obs = U( 1:Hrank , : ) * sqrt( S(1:Hrank,1:Hrank) );
% Cont = sqrt( S(1:Hrank,1:Hrank) ) * V( : , 1:Hrank)';
Obs = U( : , 1:Hrank ) * sqrt( S(1:Hrank,1:Hrank) );
Cont = sqrt( S(1:Hrank,1:Hrank) ) * V( : , 1:Hrank)';

% extract state space matrices
C = Obs(1,:);
B = Cont(:,1);
A = pinv( Obs(1:end-1,:) ) * Obs(2:end,:);

% simulate linear system model with validation inputs
x0 = zeros( size(A,1) , 1 );
y0 = 0;
x(1,:) = x0';
y(1,:) = y0';
for i = 2 : length( input_t2 )
    x(i,:) = ( A * x(i-1,:)' + B * input_u2(i-1) )';
    y(i) = C * x(i,:)';
end

%% Construct Luenberger observer on the state

% observer gain matrix

% tinder = ( 1:length(x0)/2 ) ./ ( length(x0)/2 + 1 );
% poles_desired = 0.99 .* complex( cos(tinder * 0.8*pi) , sin( tinder * 0.8*pi ) );
% poles_desired = [ poles_desired , 0.99 .* complex( cos(tinder * 0.8*pi) , -sin( tinder * 0.8*pi ) ) ];
% L = place(A',C', poles_desired ).';

% L = 1e-3 * ones( size(x0) );   % just try it out
% L = 1e-1 * pinv(C);

% sys = ss(A,B,C,[],h);
% Qn = 1e1 * eye( length(x0) ); % covariance of process noise
% Rn = 0 * eye( 1 );  % covariance of measurement noise 
% Nn = zeros( length(x0) , length(y0) );
% [~,L,~] = kalman(sys,Qn,Rn);

[num,den] = ss2tf(A,B,C,0);
[Ap,Bp,Cp,Dp] = tf2ss(num,den);
Ao = Ap'; Bo = Cp' ; Co = Bp' ; Do = Dp;
L = 0.5e-1 * Ao(:,1); % "dead-beat" controller

% L = 6e-2 * ones( size(x0) );    % Goddamnit this works the best

% set initial conditions for observer
xhat0 = x0;
% xhat0 = 1e-1 * rand(size(x0)); % initialize observer state away from origin
xhat(1,:) = xhat0';
yhat0 = Co * xhat0;
yhat(1) = yhat0';

% % simulate observer
% for i = 2 : length( input_t2 )
%     xhat(i,:) = ( A * xhat(i-1,:)' + B * input_u2(i-1) + L * ( yout2(i-1)' - yhat(i-1)' ) )';
%     yhat(i) = C * xhat(i,:)';
% end

% simulate observer (with observer canonical model)
for i = 2 : length( input_t2 )
    xhat(i,:) = ( Ao * xhat(i-1,:)' + Bo * input_u2(i-1) + L * ( yout2(i-1) - yhat(i-1) ) )';
    yhat(i) = Co * xhat(i,:)';
end

% plot comparison between observer output, real output, and model output
figure;
hold on;
plot( input_t2 , yout2( 1 , : )' , 'LineWidth' , 2);
plot( input_t2 , y , 'LineWidth' , 2);
plot( input_t2 , yhat , 'LineWidth' , 2);
hold off;

%% Construct Kalman filter observer

xhat0 = x0;
xhat(1,:) = xhat0';
yhat0 = C * xhat0;
yhat(1) = yhat0';

% initialize covariance matrices
Q = 1e2 * eye( length(x0) ); % covariance of process noise (set really high because I want filter to be aggressive)
% Q = (1:10)' * (1:10);
R = 1e1 * eye( 1 );  % covariance of measurement noise (small)
P0 = 0 * ones( size(A) );


% simulate observer
for i = 2 : length( input_t2 )
    if i == 2
        xhatk = xhat0;
        Pk = P0;
    else
        xhatk = xhatkp1;
        Pk = Pkp1;
    end
    
    % Prediction
    xkp1_ = A * xhatk + B * input_u2(i-1);
    Pkp1_ = A * Pk * A' + Q;
    
    % Measurement
    mu = 0; % mean of measurement noise
    sigma = 0.000;  % standard deviation of measurement noise
    ykp1 = yout2(i)' + normrnd(mu,sigma);   % has no noise added in

    % Correction
    Kkp1 = Pkp1_ * C' * inv( C * Pkp1_ * C' + R );
    xhatkp1 = xkp1_ + Kkp1 * ( ykp1 - C * xkp1_ );
    Pkp1 = ( eye(size(Pk) ) - Kkp1 * C ) * Pkp1_;
    
    % write state estimate to array
    xhat(i,:) = xhatkp1';
    yhat(i) = C * xhat(i,:)';

end

% plot comparison between observer output, real output, and model output
figure;
hold on;
plot( input_t2 , yout2( 1 , : )' , 'LineWidth' , 2);
plot( input_t2 , y , 'LineWidth' , 2);
plot( input_t2 , yhat , 'LineWidth' , 2);
hold off;

%% try simulating from a nonzero initial condition, using observer estimate

% try starting simulation 500 steps in
stind = 20;
input_u2_2 = input_u2( stind : end );
input_t2_2 = input_t2( stind : end );
x0_2 = xhat( stind , : )';
y0_2 = y( stind )';
x_2(1,:) = x0_2';
y_2(1) = y0_2';
for i = 2 : length( input_t2_2 )
    x_2(i,:) = ( A * x_2(i-1,:)' + B * input_u2_2( i-1) )';
    y_2(i) = C * x_2(i,:)';
end











