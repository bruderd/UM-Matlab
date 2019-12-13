%% Find T, A, B using PCA

% rename some workspace variables (for convenience)
N = sysid.params.N;
m = sysid.params.m;
n = sysid.params.n;

% take difference between before/after snapshots
% dif = sysid.koopData.Py - sysid.koopData.Px;
dif = sysid.koopData.Px;

% do pca on this difference
[ coeffs , ~ , ~ , ~ , explained , ~ ] = pca( dif );
difmean = mean( dif );  % mean value for centering data
T = coeffs';
% T = coeffs(:,1:n)';    % coordinate transformation
Tinv = pinv(T);

% transform the data (v1)
data.Px = ( sysid.koopData.Px - difmean ) * T';
data.Py = ( sysid.koopData.Py - difmean ) * T';
data.u = sysid.koopData.u;

% % transform the data usind the "double lift" (v2)
% data.x = ( sysid.koopData.Px - difmean ) * T';
% data.y = ( sysid.koopData.Py - difmean ) * T';
% data.u = sysid.koopData.u;
% % do the "double lift"
% for i = 1:length(data.x)
%     psix = sysid.lift.full( data.x(i,:)' )';
%     psiy = sysid.lift.full( data.y(i,:)' )';
%     data.Px(i,:) = psix;
%     data.Py(i,:) = psiy;     % exclude u from Py (could also use same u as Px)
% end


% learn best A, B matrices on the transformed data
cd('..');
Kvec = sysid.solve_KoopmanQP( [ data.Px , data.u ] , [ data.Py , data.u ] , sysid.lasso );
cd('findingT');
Kmtx = reshape(Kvec, [N+m,N+m]);
KT = Kmtx';    % transpose of koopman operator
A = KT( 1 : N , 1 : N );
B = KT( 1 : N , N+1 : end );


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
Tzreal = ( zreal - difmean ) * T';

% simulate the system model (v1)
% set initial condition
valdata_wzeta = sysid.get_zeta( valdata );
cd('findingT');
zeta0 = valdata_wzeta.zeta(1,:)';    % initial state with delays
z0 = sysid.lift.full( zeta0 );    % initial lifted state
Tz0 = T * ( z0 - difmean' );   % transformed initial state

% simulate lifted linear model
tsim = treal;
usim = ureal;
ysim = zeros( size( yreal ) ); % preallocate
Tzsim = zeros( size( zreal ) ); % preallocate
ysim(1,:) = yreal(1,:); % initialize
Tzsim(1,:) = Tz0';        % initialize
for j = 1 : length(treal)-1
    Tzsim(j+1,:) = ( A * Tzsim(j,:)' + B * usim(j,:)' )';
    ysim(j+1,:) = ( sysid.model.C * ( Tinv * Tzsim(j+1,:)' + difmean' ) )';
end

% % simulate the system model (v2)
% % set initial condition
% valdata_wzeta = sysid.get_zeta( valdata );
% cd('findingT');
% zeta0 = valdata_wzeta.zeta(1,:)';    % initial state with delays
% z0 = sysid.lift.full( zeta0 );    % initial lifted state
% Tz0 = T * ( z0 - difmean' );   % transformed initial state
% Z0 = sysid.lift.full( Tz0 );   % lift the transformed state
% 
% % simulate lifted linear model
% tsim = treal;
% usim = ureal;
% ysim = zeros( size( yreal ) ); % preallocate
% Zsim = zeros( size( zreal ) ); % preallocate
% ysim(1,:) = yreal(1,:); % initialize
% Zsim(1,:) = Z0';        % initialize
% for j = 1 : length(treal)-1
%     Zsim(j+1,:) = ( A * Zsim(j,:)' + B * usim(j,:)' )';
%     ysim(j+1,:) = ( sysid.model.C * ( Tinv * sysid.model.C * Zsim(j+1,:)' + difmean' ) )';
% end


