%% find T, A, and B

% rename some workspace variables
N = sysid.params.N;
m = sysid.params.m;

% Use linear A and B matrices for initial condition
T0 = eye( N );
A0 = sysid.model.A; % eye(N);
B0 = sysid.model.B; % zeros( N,m );

% Vectorize and stack matricies for initial decision variable
arg0 = [ T0(:) ; A0(:) ; B0(:) ];


% solve optimization problem
opts = optimoptions( 'fmincon' , 'Algorithm' , 'interior-point' , 'Display' , 'iter' , 'MaxFunctionEvaluations' , 1e4 );
argmin = fmincon( @(arg) costfun( arg , sysid.koopData , sysid.params ) , arg0 , [] , [] , [] , [] , [] , [] , @(arg) nonlcon( arg , sysid.params ) , opts );

% convert solution back into matrices
T = reshape( argmin( 1 : N^2 ) , [N,N] );
Tinv = inv(T);  % inverse transformation
A = reshape( argmin( N^2 + 1 : 2*N^2 ) , [N,N] );
B = reshape( argmin( 2*N^2 + 1 : end ) , [N,m] );

%% simulate new system against validation data

% rename validation data struct
valdata = sysid.valdata{1};

% shift real data so delays can be accounted for
index0 = sysid.params.nd + 1;  % index of the first state
treal = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
yreal = valdata.y(index0 : end , :);
ureal = valdata.u(index0 : end , :);
cd('..');   % change directories so you can call sysid.get_zeta function
[ ~ , zetareal ] = sysid.get_zeta( valdata );
zreal = zeros( size( zetareal , 2 ) , sysid.params.N );
for i = 1 : size( zetareal , 1 )
    zreal(i,:) = sysid.lift.full( zetareal(i,:)' );
end


% set initial condition
valdata_wzeta = sysid.get_zeta( valdata );
cd('findingT');
zeta0 = valdata_wzeta.zeta(1,:)';    % initial state with delays
z0 = sysid.lift.full( zeta0 );    % initial lifted state
Tz0 = T * z0;   % transformed initial state

% simulate lifted linear model
tsim = treal;
usim = ureal;
ysim = zeros( size( yreal ) ); % preallocate
zsim = zeros( size( zreal ) ); % preallocate
ysim(1,:) = yreal(1,:); % initialize
Tzsim(1,:) = Tz0';        % initialize
for j = 1 : length(treal)-1
    Tzsim(j+1,:) = ( A * Tzsim(j,:)' + B * usim(j,:)' )';
    ysim(j+1,:) = ( sysid.model.C * Tinv * Tzsim(j+1,:)' )';
end
