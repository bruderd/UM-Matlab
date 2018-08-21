function [ sysid_data, val_data ] = partition_trimdata(val_len)
%partition_data: Partitions a trimmed data file into sysid and validation
%components. Will take the first 'val_len' seconds of data for the
%validation component, and the rest for the sysid component
%   INPUTS
%       val_len - a constant defining the length of the validation data in
%       seconds
%

% Prompt user to identify data file
disp(['Please select .mat file corresponding to trimmed data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
rawdata = load([data_path, data_file]);


% find index of first data point past val_len
for index = 1:length(rawdata.t)
    if rawdata.t(index) >= val_len
        break;
    end
end

%% val data
val_data = struct;
val_data.t = rawdata.t(1:index,:);
val_data.x = rawdata.x(1:index,:);
val_data.u = rawdata.u(1:index,:);

% save validation data file
t = val_data.t;
x = val_data.x;
u = val_data.u;
val_fname = ['trimdataFiles', filesep, 'validation', filesep, data_name, '_val_', num2str(val_len), 's.mat'];
save(val_fname, 't', 'x', 'u');

%% sysid data
sysid_data = struct;
sysid_data.t = rawdata.t(index+1 : end) - rawdata.t(index+1);
sysid_data.x = rawdata.x(index+1 : end, :);
sysid_data.u = rawdata.u(index+1 : end, :);

% save sysid data file
t = sysid_data.t;
x = sysid_data.x;
u = sysid_data.u;
sysid_fname = ['trimdataFiles', filesep, 'sysid', filesep, data_name, '_sysid_', num2str(val_len), 's.mat'];
save(sysid_fname, 't', 'x', 'u');


end

