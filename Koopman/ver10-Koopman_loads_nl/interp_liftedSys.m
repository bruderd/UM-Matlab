function modb = interp_liftedSys(a,b,c)
% interp_liftedSys
%
% Finds an approximate Koopman operator corresponding to a system with a
% load attached whose mass is in between the masses of two other Koopman
% models.

%% Load in data and model files

% Prompt user to identify data file
[data_file,data_path] = uigetfile( 'dataFiles' );
matcontents = load([data_path, data_file]); % must be a .mat file
data = matcontents.data;

% Prompt user to identify model files
[moda_file,moda_path] = uigetfile( ['models', filesep] , ['Select model for load = ' , num2str(a)] );
[modc_file,modc_path] = uigetfile( ['models', filesep] , ['Select model for load = ' , num2str(c)] );
moda = load([moda_path , moda_file]);
modc = load([modc_path , modc_file]);

% Get the lifting function (must be the same for both models)
cd( ['liftingFunctions'] );
lift = str2func( [ 'lift_' , moda.params.systemName ] );
cd(['..']);

%% Build matrices

theta = ( b-a ) / ( c-a );
[x,y,u] = deal(data.snapshotPairs.zeta_x, data.snapshotPairs.zeta_y, data.snapshotPairs.u);

for i = 1:length( x(:,1) )
    psix = lift( x(i,:)' )';
%     psiy = lift( y(i,:)' )';
    Px(i,:) = [ psix , u(i,:) ];
%     Py(i,:) = [ psiy , zeros(1,p) ];     % exclude u from Py (could also use same u as Px
end

Omega = (1-theta) * Px * moda.U + theta * Px * modc.U;

%% Solve for new Koopman operator and identify new linear model
U = pinv( Px ) * Omega;

UT = U';    % transpose of koopman operator

A = UT( 1 : moda.params.N , 1 : moda.params.N );
B = UT( 1 : moda.params.N , moda.params.N+1 : end );

% if strcmp( params.basisID, 'poly' )
%     C = [zeros(params.ny , 1) , eye(params.ny) , zeros( params.ny , params.N - params.ny - 1 )];   % if poly we want to skip over first element of lifted state which is "1"
% else
    C = [eye(moda.params.ny), zeros(moda.params.ny , moda.params.N - moda.params.ny)];   % C selects the first ny entries of the lifted state (so output can be different than state)
% end

% matrix to recover the state with delays, i.e. zeta
Cz = [eye(moda.params.nzeta), zeros(moda.params.nzeta , moda.params.N - moda.params.nzeta)];

%% Define output struct
modb = struct;

modb.A = A;  
modb.B = B;  
modb.C = C;
modb.Cz = Cz;
modb.sys = ss(A,B,C,0, moda.params.Ts);  % discrete state space system object
modb.params = moda.params;    % save system parameters as part of system struct
modb.U = U;  % koopman operator matrix

