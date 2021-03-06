% sysid_unl_main
%
% Creates a sysid_unl class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)


%% gather training data (need to prepare data file before running this)

% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );


%% construct sysid class
options = {};
sysid_unl = sysid_unl( data4sysid, options ,...
    'obs_type' , { 'poly' } ,...
    'obs_degree' , [ 2 ] ,...
    'snapshots' , Inf ,...
    'lasso' , [ Inf ] ,...
    'delays' , 0 ,...
    'isupdate' , false,...
    'armclass' , [] ,...
    'liftinput' , 1 );

%% validate model(s)
% could also manually do this for one model at a time

results = cell( size(sysid_unl.candidates) );    % store results in a cell array
err = cell( size(sysid_unl.candidates) );    % store error in a cell array 

if iscell(sysid_unl.candidates)
    for i = 1 : length(sysid_unl.candidates)
        [ results{i} , err{i} ] = sysid_unl.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = sysid_unl.valNplot_model;
end

%% calculate the  no-input-model error on the training data

% [ e , sysid_unl.snapshotPairs ] = sysid_unl.get_e( sysid_unl.snapshotPairs );   % use snapshot pairs
[ e , sysid_unl.traindata ] = sysid_unl.get_e( sysid_unl.traindata );   % use timeseries data

% get a reduced order version of nu
[ nu , Beta , sysid_unl.traindata ] = sysid_unl.get_nu( sysid_unl.traindata );

% also calculate nu for the validation data
for i = 1 : size( sysid_unl.valdata , 2 )
    [ ~ , sysid_unl.valdata{i} ] = sysid_unl.get_e( sysid_unl.valdata{i} );
    sysid_unl.valdata{i}.nu = sysid_unl.valdata{i}.e * Beta';
end

%% train neural network to map from nu to u

noise = 0.1 * ( rand( size(sysid_unl.traindata.nu) ) - 0.5 );   % add noise to make the neural network more robust to disturbances
[ sysid_unl.e2u.nnet , sysid_unl.e2u.fun ] = train_nnet( [sysid_unl.traindata.nu + noise, sysid_unl.traindata.y(1:end-1,:) ] , sysid_unl.traindata.u( 1 : size( nu , 1 ) , : )  );           %train from nu and state
% [ sysid_unl.e2u.nnet , sysid_unl.e2u.fun ] = train_nnet( [sysid_unl.traindata.nu ] , sysid_unl.traindata.u( 1 : size( nu , 1 ) , :) );    % train from nu
% [ sysid_unl.e2u.nnet , sysid_unl.e2u.fun ] = train_nnet( [sysid_unl.traindata.e , sysid_unl.traindata.y(1:end-1,:) ] , sysid_unl.traindata.u( 1 : size( e , 1 ) , : )  );           %train from e and state

%% modify model so that it will work with mpc

sysid_unl.model.Beta = Beta;
sysid_unl.model.Beta_pinv = pinv(Beta);
sysid_unl.model.params.mnu = size( nu , 2 );
sysid_unl.params.mnu = size( nu , 2 );
sysid_unl.params.nzx = size( sysid_unl.basis.zx , 1 );
sysid_unl.params.nzu = size( sysid_unl.basis.zu , 1 );

sysid_unl.params.NLinput = sysid_unl.e2u.fun;

% %% train neural network to map from e to u
% 
% [ sysid_unl.e2u.nnet , sysid_unl.e2u.fun ] = train_nnet( sysid_unl.traindata.e , sysid_unl.traindata.u(1:end-(1+sysid_unl.params.nd),:)  );
% 
% %% modify model so that it will work with mpc
% 
% sysid_unl.model.Beta = speye( sysid_unl.model.params.N );
% sysid_unl.model.B = speye( sysid_unl.model.params.N );
% % sysid_unl.model.params.m = sysid_unl.model.params.N;
% % sysid_unl.params.m = sysid_unl.params.N;
% 
% sysid_unl.params.NLinput = sysid_unl.e2u.fun;

%NEED TO MODIFY BELOW THIS LINE--------------------------------------------      


%% save model(s)

% You do this based on the validation results.
% Call this function:
%   sysid.save_class( model_id )

