function sliceNdice_2sVals(val_len , num_vals)
%sliceNdice_2sVals: Creates num_vals number of validation trials of length
%val_len by taking random chunks of a longer trimmed data file
%   INPUTS
%       val_len - a constant defining the length of the validation data in
%       seconds
%       num_vals - number of slices to take
%
%   NOTE: This only differs from 'sliceNdice_trimmed' in that it saves the
%   files in the folder 'waves_2sVals'

% Prompt user to identify data file
disp(['Please select .mat file corresponding to trimmed data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
rawdata = load([data_path, data_file]);


% choose random starting points for all the validation trials
start_time = ceil( ( rawdata.t(end) - val_len ) * rand( num_vals , 1 ) );
start_ind = zeros( num_vals , 1 );
for i = 1 : num_vals
    for j = 1:length(rawdata.t)
        if rawdata.t(j) >= start_time(i)
            start_ind(i) = j;
            break;
        end
    end
end

% find the number of steps in val_len
for index = 1:length(rawdata.t)
    if rawdata.t(index) >= val_len
        break;
    end
end

%% save the data for each trial
val_data = struct;
for i = 1 : num_vals
    val_data.t = rawdata.t( start_ind(i) : start_ind(i)+index , : ) - rawdata.t( start_ind(i) );
    val_data.x = rawdata.x( start_ind(i) : start_ind(i)+index , : );
    val_data.u = rawdata.u( start_ind(i) : start_ind(i)+index , : );
    
    % save validation data file
    t = val_data.t;
    x = val_data.x;
    u = val_data.u;
    val_fname = ['trimdataFiles', filesep, 'waves_2sVals', filesep, data_name, '_val_', num2str(val_len), 's_' , num2str(i) , '.mat'];
    save(val_fname, 't', 'x', 'u');
end