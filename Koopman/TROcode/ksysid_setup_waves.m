% ksysid_setup_waves
%
% Creates a sysid class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)
%
% Modified to run on waves with no GUI


%% gather training data (need to prepare data file before running this)

% load in data file(s)
data4sysid = load( [ 'datafiles' , filesep , 'RSS2019-robot_train-16_val-3_2020-05-06_18-42.mat' ] );   % older robot data

%% construct sysid class
ksysid = ksysid( data4sysid, ...
        'model_type' , 'linear' ,...    % model type (linear or nonlinear)
        'obs_type' , { 'poly' } ,...    % type of basis functions
        'obs_degree' , [ 4 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ 0:0.5:10 ] ,...           % L1 regularization term
        'delays' , 1 );                 % Number of state/input delays

    
%% train model(s)
ksysid = ksysid.train_models;


%% validate model(s) (Skip since no GUI in waves)
% could also manually do this for one model at a time

% results = cell( size(ksysid.candidates) );    % store results in a cell array
% err = cell( size(ksysid.candidates) );    % store error in a cell array 
% 
% if iscell(ksysid.candidates)
%     for i = 1 : length(ksysid.candidates)
%         [ results{i} , err{i} ] = ksysid.valNplot_model( i );
%     end
% else
%     [ results{1} , err{1} ] = ksysid.valNplot_model;
% end
    

%% save model(s)

% You do this based on the validation results.
% Call this function:
  ksysid.save_class( [] , [ 'systems' , filesep , 'waves' ] );
