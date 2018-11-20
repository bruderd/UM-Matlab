function [t,x,u] = trim_data(tstart, tend)
%trim_data: trim away the beginning and end of raw data file
%   INPUTS
%       tstart - constant specifying time at which trimmed data should
%       start (but will be the new zero)
%       tend - constant specifying time at which the trimmed data should
%       end

% Prompt user to identify data file
disp(['Please select .mat file corresponding to raw data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
rawdata = load([data_path, data_file]);


% find index of first data point past tstart
for i_tstart = 2:length(rawdata.t)  % first point can be shitty so skip it
    if rawdata.t(i_tstart) >= tstart
        break;
    end
end

% find index of first data point past tend
for i_tend = 1:length(rawdata.t)
    if rawdata.t(i_tend) >= tend
        break;
    end
end

% define trimmed data series
t = rawdata.t(i_tstart : i_tend) - rawdata.t(i_tstart);
x = rawdata.x(i_tstart : i_tend , :);
u = rawdata.u(i_tstart : i_tend , :);

% save trimmed data files
% val_fname = ['trimdataFiles', filesep, data_name, '_200s.mat'];
val_fname = ['trimdataFiles', filesep, data_name, '_trim.mat'];
save(val_fname, 't', 'x', 'u');

end
