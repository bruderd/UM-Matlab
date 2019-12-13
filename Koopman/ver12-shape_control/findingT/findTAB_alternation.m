%% find T, A, and B

% USER EDIT
niter = 6; % number of iterations/alternations

% rename some workspace variables (for convenience)
N = sysid.params.N;
m = sysid.params.m;

% initialize cell arrays to contain A,B,T at each iteration
Aall = cell( niter , 1 );
Ball = cell( niter , 1 );
Tall = cell( niter , 1 );

% Use linear A and B matrices for initial condition
T0 = eye( N );  % diag( rand(N,1) ); 
T0(1,:) = [ -3 0 0 2 0 0 ];
A0 = sysid.model.A; % eye(N);
B0 = sysid.model.B; % zeros( N,m );

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
%     argmin = fmincon( @(arg) costfun_Tonly( arg , Aall{i-1} , Ball{i-1} , data_i , sysid ) , arg0 , [] , [] , [] , [] , [] , [] , @(arg) nonlcon( arg , sysid.params ) , opts );
    argmin = fmincon( @(arg) costfun_Tonly( arg , Aall{i-1} , Ball{i-1} , sysid.koopData , sysid) , arg0 , [] , [] , [] , [] , [] , [] , @(arg) nonlcon( arg , sysid.params ) , opts );

    % convert solution back into matrix
    Tall{i} = reshape( argmin( 1 : N^2 ) , [N,N] );
%     T_trans = T0';
%     for j = 2 : i
%         T_trans = T_trans * Tall{j}';
%     end
    T_trans = Tall{i}';
    
    % transform data 
    data_i.Px = sysid.koopData.Px * T_trans;
    data_i.Py = sysid.koopData.Py * T_trans;
    data_i.u = sysid.koopData.u;
    
    % Solve for new A,B matrices
    cd('..');
    Kvec = sysid.solve_KoopmanQP( [ data_i.Px , data_i.u ] , [ data_i.Py , data_i.u ] , sysid.lasso );
    cd('findingT');
    Kmtx = reshape(Kvec, [N+m,N+m]);
    KT = Kmtx';    % transpose of koopman operator
    Aall{i} = KT( 1 : N , 1 : N );
    Ball{i} = KT( 1 : N , N+1 : end );

end

% pull out final matrices
T = T_trans';   % want all the transformations put together
Tinv = pinv(T);
A = Aall{end};
B = Ball{end};

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
Tzsim = zeros( size( zreal ) ); % preallocate
ysim(1,:) = yreal(1,:); % initialize
Tzsim(1,:) = Tz0';        % initialize
for j = 1 : length(treal)-1
    Tzsim(j+1,:) = ( A * Tzsim(j,:)' + B * usim(j,:)' )';
    ysim(j+1,:) = ( sysid.model.C * Tinv * Tzsim(j+1,:)' )';
end