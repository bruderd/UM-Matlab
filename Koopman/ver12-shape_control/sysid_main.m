% sysid_main
%
% Creates a sysid class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)
% HEAVILY MODIFIED 8/17/2019. Still have commented out peices but should
% delete most of the code later...


%% gather training data (need to prepare data file before running this)

% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );


%% construct sysid class
sysid = sysid( data4sysid, ...
        'obs_type' , { 'poly' } ,...
        'obs_degree' , [ 2 ] ,...
        'snapshots' , Inf ,...
        'lasso' , [ Inf ] ,...
        'delays' , 0 ,...
        'isupdate' , false,...
        'armclass' , [] ); % this is only needed if you want armshape observables

    
%% train model(s)
sysid = sysid.train_models;


%% validate model(s)
% could also manually do this for one model at a time

results = cell( size(sysid.candidates) );    % store results in a cell array
err = cell( size(sysid.candidates) );    % store error in a cell array 

if iscell(sysid.candidates)
    for i = 1 : length(sysid.candidates)
        [ results{i} , err{i} ] = sysid.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = sysid.valNplot_model;
end
    

%% save model(s)

% You do this based on the validation results.
% Call this function:
%   sysid.save_class( model_id )


