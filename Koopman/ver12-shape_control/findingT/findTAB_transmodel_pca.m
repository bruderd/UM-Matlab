%% find transformed model (i.e. "double lift" using pca)

% USER EDIT
niter = 8; % number of iterations/alternations

% rename some workspace variables (for convenience)
N = sysid.params.N;
m = sysid.params.m;
n = sysid.params.n;

% take difference between before/after snapshots
dif = sysid.koopData.Py - sysid.koopData.Px;
% dif = sysid.koopData.Px;

% do pca on this difference
[ coeffs , ~ , ~ , ~ , explained , ~ ] = pca( dif );
difmean = mean( dif );  % mean value for centering data
% T = coeffs';
T = coeffs(:,1:n)';    % coordinate transformation
Tinv = pinv(T);


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