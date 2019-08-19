% sysid_main
%
% Creates a sysid class and walks through all of the steps of building a
% model from data, validating its performance, and saving it (if desired)
% HEAVILY MODIFIED 8/17/2019. Still have commented out peices but should
% delete most of the code later...

% (USER EDIT HERE)
isupdate = false; % put c

%% gather training data (need to prepare data file before running this)

% load in data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );

% % separate into training and validation data
% data4train = data4sysid.train;
% data4val = data4sysid.val;
% 
% % create sysid class from data
% sysid = sysid( data4train{1} , data4train{1}.params );
% 
% % merge the training data into a single big file
% data4train_merged = sysid.merge_trials( data4train );
% 
% % scale data to be in range [-1 , 1]
% [ traindata , sysid ] = sysid.get_scale( data4train_merged );
% valdata = cell( size( data4val ) );  
% for i = 1 : length( data4val )
%     valdata{i} = sysid.scale_data( data4val{i} );
% end

%% construct sysid class
sysid = sysid( data4sysid, ...
        'obs_type' , { 'poly' } ,...
        'obs_degree' , [ 2 ] ,...
        'snapshots' , Inf ,...
        'lasso' , [ 1 , 2 ] ,...
        'delays' , 1 );

%% define a set of observables

% % (USER EDIT HERE)
% type = { 'armshape' , 'poly' }; % types of functions to include
% degree = [ 1 , 3 ]; % degree/complexity of each type of function
% 
% sysid = sysid.def_observables( type , degree );


%% learn a Koopman model

% % (USER EDIT HERE)
% numPairs = Inf;  % number of snapshot pairs to be taken (Inf will take maximum amount)
% lasso = []; % lasso L1 regularization parameter. ( empty array for default value, which gives least-squares solution basically )
% 
% snapshotPairs = sysid.get_snapshotPairs( traindata , numPairs );
[ koopData , K ] = sysid.get_Koopman( sysid.snapshotPairs , sysid.lasso );
[ model , sysid ] = sysid.get_model( koopData );

%% validate model

results = cell( size(sysid.valdata) );    % store results in a cell array
err = cell( size(sysid.valdata) );    % store error in a cell array 
for i = 1 : length(sysid.valdata)
    results{i} = sysid.val_model( model , sysid.valdata{i} );
    err{i} = sysid.get_error( results{i}.sim , results{i}.real );
    sysid.plot_comparison( results{i}.sim , results{i}.real );
end

%% save (or don't save) sysid class, model, and training data

saveModel = questdlg('Would you like to save this model?');
if strcmp( saveModel , 'Yes' )
    model = sysid.save_class( model , isupdate);
end




