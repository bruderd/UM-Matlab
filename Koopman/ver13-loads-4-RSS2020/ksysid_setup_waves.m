% ksysid_setup_waves
%
% Creates a sysid class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)
%
% Modified to work on waves by removing things that require a gui such as
% 'uigetfile', ...


%% gather training data (need to prepare data file before running this)

% % load in data file(s)
% data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_diffvaldata_train-5_val-5_2020-01-03_12-03.mat' ]);  % Xun's first loaded tests
% data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_7markers_0p5s-ramp_train-6_val-6_2020-01-06_14-16.mat' ]);  % Xun's loaded tests with faster ramp (0.5s)
% data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_3markers_0p5-1s-ramps_train-11_val-11_2020-01-06_18-44.mat' ]);  % Xun's loaded tests with slow and faster ramp (0.5s), only 3 rigid bodies for state
% data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_1markers_0p5-1s-ramps_train-11_val-11_2020-01-07_12-45.mat' ]);  % same data as above, but just the end effector
% data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_3marks_0p5-1s-ramps_train-22_val-22_2020-01-09_13-25.mat' ]);    % includes more loads from 0-500g for the 0.5 and 1 sec ramps
data4sysid = load([ 'datafiles' , filesep , 'softarm_3mods_3marks_0p5-1s-ramps_trim_train-22_val-22_2020-01-10_20-11.mat' ]);   % same as last file, but with start/end cut off from each trial so avoid weird bad stuff

%% construct sysid class
ksysid = ksysid( data4sysid, ...
        'model_type' , 'linear' ,...    % model type (linear or nonlinear)
        'obs_type' , { 'poly' } ,...    % type of basis functions
        'obs_degree' , [ 2 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ Inf ] ,...           % L1 regularization term
        'delays' , 1 ,...               % Numer of state/input delays
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
%     

%% save model(s)

ksysid.save_class;


