function modb = interp_multi_liftedSys( tq , t )
% interp_liftedSys
%
% Finds an approximate Koopman operator corresponding to a system with a
% load attached whose mass is in between the masses of two other Koopman
% models.
%
% tq = load for the unknown model
% t = vector of loads for each of the unknown models

%% Load in data and model files

% Prompt user to identify data file
[data_file,data_path] = uigetfile( 'dataFiles' , 'Select a date file which will provide the sample points (something near the desired model)' );
matcontents = load([data_path, data_file]); % must be a .mat file
data = matcontents.data;

% Prompt user to identify model files
for i = 1 : length(t)
    [mod_file,mod_path] = uigetfile( ['models', filesep] , ['Select model for load = ' , num2str(t(i))] );
    mod{i} = load([mod_path , mod_file]);
end

% Get the lifting function (must be the same for all models)
cd( ['liftingFunctions'] );
lift = str2func( [ 'lift_' , mod{1}.params.systemName ] );
cd( ['..'] );

%% Build matrices

% snapshot pairs. We will only use the x's
[x,y,u] = deal(data.snapshotPairs.zeta_x, data.snapshotPairs.zeta_y, data.snapshotPairs.u);

Px = zeros( length(x(:,1)) , mod{1}.params.Np );
for i = 1:length( x(:,1) )
    psix = lift( x(i,:)' )';
%     psiy = lift( y(i,:)' )';
    Px(i,:) = [ psix , u(i,:) ];
%     Py(i,:) = [ psiy , zeros(1,p) ];     % exclude u from Py (could also use same u as Px
end

% Interpolation (cubic)
Py = zeros( length(x(:,1)) , mod{1}.params.Np );
for i = 1 : length( x(:,1) )
    for j = 1 : length(t)
        psiq(j,:) = Px(i,:) * mod{j}.U ;
    end
    Py(i,:) = interp1( t , psiq , tq , 'linear');
end

%% Solve for new Koopman operator and identify new linear model
U = pinv( Px ) * Py;

UT = U';    % transpose of koopman operator

A = UT( 1 : mod{1}.params.N , 1 : mod{1}.params.N );
B = UT( 1 : mod{1}.params.N , mod{1}.params.N+1 : end );

% if strcmp( params.basisID, 'poly' )
%     C = [zeros(params.ny , 1) , eye(params.ny) , zeros( params.ny , params.N - params.ny - 1 )];   % if poly we want to skip over first element of lifted state which is "1"
% else
    C = [eye(mod{1}.params.ny), zeros(mod{1}.params.ny , mod{1}.params.N - mod{1}.params.ny)];   % C selects the first ny entries of the lifted state (so output can be different than state)
% end

% matrix to recover the state with delays, i.e. zeta
Cz = [eye(mod{1}.params.nzeta), zeros(mod{1}.params.nzeta , mod{1}.params.N - mod{1}.params.nzeta)];

%% Define output struct
modb = struct;

modb.A = A;  
modb.B = B;  
modb.C = C;
modb.Cz = Cz;
modb.sys = ss(A,B,C,0, mod{1}.params.Ts);  % discrete state space system object
modb.params = mod{1}.params;    % save system parameters as part of system struct
modb.U = U;  % koopman operator matrix