% ksysid_setup
%
% Creates a sysid class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)


%% gather training data (need to prepare data file before running this)

% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );


%% construct sysid class
ksysid = ksysid( data4sysid, ...
        'model_type' , 'linear' ,...    % model type (linear or nonlinear)
        'obs_type' , { 'poly' } ,...    % type of basis functions
        'obs_degree' , [ 2 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ Inf ] ,...           % L1 regularization term
        'delays' , 0 ,...               % Numer of state/input delays
        'loaded' , true);             % Does system include loads?
    
disp(['Number of basis functions: ' , num2str( 2 * ksysid.params.N ) ]);
    
%% basis dimensional reduction (beta)

disp('Performing dimensional reduction...');
Px = ksysid.lift_snapshots( ksysid.snapshotPairs );
ksysid = ksysid.get_econ_observables( Px );
disp(['Number of basis functions: ' , num2str( 2 * ksysid.params.N ) ]);
clear Px;
    
%% train model(s)
ksysid = ksysid.train_models;


%% validate model(s)
% could also manually do this for one model at a time

results = cell( size(ksysid.candidates) );    % store results in a cell array
err = cell( size(ksysid.candidates) );    % store error in a cell array 

if iscell(ksysid.candidates)
    for i = 1 : length(ksysid.candidates)
        [ results{i} , err{i} ] = ksysid.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = ksysid.valNplot_model;
end
    

%% save model(s)

% You do this based on the validation results.
% Call this function:
%   ksysid.save_class( )


