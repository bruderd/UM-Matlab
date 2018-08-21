function [ sysid_data, val_data ] = partition_data
%partition_data: Partitions a trimmed data file into sysid and validation
%components

% Prompt user to identify data file
disp(['Please select .mat file corresponding to raw data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
rawdata = load([data_path, data_file]);



% Use 7/10 of data for sysid, 1/10 for validation, don't use beginning or
% ending 1/10

len = length(rawdata.t);
frac = floor(len/10);

%% sysid data
sysid_data = struct;
ssid = frac+1;  % starting index
esid = 8*frac;  % ending index
isid = ssid : esid; % range used
sysid_data.t = rawdata.t(isid,:) - rawdata.t(ssid);
sysid_data.x = rawdata.x(isid,:);
sysid_data.u = rawdata.u(isid,:);

% save sysid data file
t = sysid_data.t;
x = sysid_data.x;
u = sysid_data.u;
sysid_fname = ['rawdataFiles', filesep, data_name, '_sysid.mat'];
save(sysid_fname, 't', 'x', 'u');

%% val data
val_data = struct;
sval = esid+1;  % starting index
ival = sval : 9*frac; % range used
val_data.t = rawdata.t(ival,:) - rawdata.t(sval);
val_data.x = rawdata.x(ival,:);
val_data.u = rawdata.u(ival,:);

% save validation data file
t = val_data.t;
x = val_data.x;
u = val_data.u;
val_fname = ['rawdataFiles', filesep, data_name, '_val.mat'];
save(val_fname, 't', 'x', 'u');


end

