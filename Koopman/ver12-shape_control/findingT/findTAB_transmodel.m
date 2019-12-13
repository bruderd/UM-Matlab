%% find transformed model (i.e. "double lift")

% USER EDIT
niter = 2; % number of iterations/alternations

% rename some workspace variables (for convenience)
N = sysid.params.N;
m = sysid.params.m;
n = sysid.params.n;

% initialize cell arrays to contain A,B,T at each iteration
Aall = cell( niter , 1 );
Ball = cell( niter , 1 );
Tall = cell( niter , 1 );

% Use linear A and B matrices for initial condition
T0 = eye( N );  % diag( rand(N,1) ); 
T0 = T0( 1 : n , : );   % make it n x N
A0 = eye(n);
B0 = zeros( n , m );

Tall{1} = T0;
Aall{1} = A0;
Ball{1} = B0;

% initializations
arg0 = T0(:);
data_i.Px = sysid.koopData.Px;
data_i.Py = sysid.koopData.Py;
data_i.u = sysid.koopData.u;

% Solve for T,A,B through alternation
for i = 2 : niter
    
    % Vectorize T for initial decision variable
    arg0 = Tall{i-1}(:);
    
    % solve optimization problem
    opts = optimoptions( 'fmincon' , 'Algorithm' , 'interior-point' , 'Display' , 'iter' ); % 'MaxFunctionEvaluations' , 1e4 );
    argmin = fmincon( @(arg) costfun_transmodel( arg , Aall{i-1} , Ball{i-1} , sysid.koopData , sysid) , arg0 , [] , [] , [] , [] , [] , [] , @(arg) nonlcon_transmodel( arg , sysid.params ) , opts );

    % convert solution back into matrix
    Tall{i} = reshape( argmin , [n,N] );
    T_trans = Tall{i}';
    
    % transform data 
    data_i.Px = sysid.koopData.Px * T_trans;
    data_i.Py = sysid.koopData.Py * T_trans;
    data_i.u = sysid.koopData.u;
    
    % Solve for new A,B matrices
    sysid.params.N = n;     % change N locally 
    cd('..');
    Kvec = sysid.solve_KoopmanQP( [ data_i.Px , data_i.u ] , [ data_i.Py , data_i.u ] , sysid.lasso );
    cd('findingT');
    Kmtx = reshape(Kvec, [n+m,n+m]);
    KT = Kmtx';    % transpose of koopman operator
    Aall{i} = KT( 1 : n , 1 : n );
    Ball{i} = KT( 1 : n , n+1 : end );
    sysid.params.N = N;     % change N back to what it was

end

% pull out final matrices
T = T_trans';   % want all the transformations put together
Tinv = pinv(T);
A = Aall{end};
B = Ball{end};

%% learn new koopman model on transformed state v = T z

% new unlifted snapshot pairs
snapshotPairs.alpha = sysid.koopData.Px * T';
snapshotPairs.beta = sysid.koopData.Py * T';
snapshotPairs.u = sysid.koopData.u;

cd('..');
[ vkoopData , vK ] = sysid.get_Koopman( snapshotPairs );
vmodel = sysid.get_model( vkoopData );
cd('findingT');

%% simulate new system against validation data

% rename validation data struct
valdata = sysid.valdata{1};

% shift real data so delays can be accounted for
index0 = sysid.params.nd + 1;  % index of the first state
treal = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
% yreal = valdata.y(index0 : end , :);
ureal = valdata.u(index0 : end , :);
cd('..');   % change directories so you can call sysid.get_zeta function
[ ~ , zetareal ] = sysid.get_zeta( valdata );
zreal = zeros( size( zetareal , 2 ) , sysid.params.N );
for i = 1 : size( zetareal , 1 )
    zreal(i,:) = sysid.lift.full( zetareal(i,:)' );
end
yreal = zreal * T'; % transform from y to v


% set initial condition
cd('findingT');
z0 = zreal(1,:)';    % initial lifted state
v0 = T * z0;   % transformed initial state
z0 = sysid.lift.full( v0 );

% simulate lifted linear model
tsim = treal;
usim = ureal;
ysim = zeros( size( yreal ) ); % preallocate
zsim = zeros( size( zreal ) ); % preallocate
ysim(1,:) = v0'; % initialize
zsim(1,:) = z0';        % initialize
for j = 1 : length(treal)-1
    zsim(j+1,:) = ( vmodel.A * zsim(j,:)' + vmodel.B * usim(j,:)' )';
    ysim(j+1,:) = ( vmodel.C * zsim(j+1,:)' )';
end




