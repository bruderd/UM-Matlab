% make_train_val_cellarrays
%
%   This is only used to prep data for compaibility with TROcode

%% Training data

% Prompt user to select all the training files
disp(['Please select .mat file(s) corresponding to training data.']);
[data_file,data_path] = uigetfile('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver8-Koopman_MPC\readinData\trimdataFiles_v2_xps' , 'MultiSelect' , 'on');

numtrain = length( data_file );
train = cell( 1 , numtrain );
for i = 1 : length( data_file )
%     [data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
    rawdata = load([data_path, data_file{i}]);
    train{i} = rawdata;  
end

%% Validation data

% Prompt user to select all the training files
disp(['Please select .mat file(s) corresponding to validation data.']);
[data_file,data_path] = uigetfile('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver8-Koopman_MPC\readinData\trimdataFiles_v2_xps' , 'MultiSelect' , 'on');

numtrain = length( data_file );
val = cell( 1 , numtrain );
for i = 1 : length( data_file )
%     [data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
    rawdata = load([data_path, data_file{i}]);
    val{i} = rawdata;  
end