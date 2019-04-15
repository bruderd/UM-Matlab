function [t,x,u] = trim_laserData(tstart, tend)
%trim_data: trim away the beginning and end of raw data file from the
%cartesion laser tracking experiments
%   INPUTS
%       tstart - constant specifying time at which trimmed data should
%       start (but will be the new zero)
%       tend - constant specifying time at which the trimmed data should
%       end

% Prompt user to identify data file
disp(['Please select .mat file corresponding to raw data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
foo = load([data_path, data_file]);
rawdata = foo.sysidData;

%% treat the data
% Remove all NaNs from the data
rawdata.T(any(isnan(rawdata.Y), 2), :) = [];
rawdata.Y(any(isnan(rawdata.Y), 2), :) = [];
rawdata.U(any(isnan(rawdata.Y), 2), :) = [];

% interpolate over outliers (using cubic spline interp) (careful, this can
% cause problems. MAY NEED TO TRIM DATA AGAIN WITHOUT THIS!!!)
% rawdata.Y = filloutliers(rawdata.Y,'linear',1);

%% trim the data
% find index of first data point past tstart
for i_tstart = 2:length(rawdata.T)  % first point can be shitty so skip it
    if rawdata.T(i_tstart) >= tstart
        break;
    end
end

% find index of first data point past tend
for i_tend = 1:length(rawdata.T)
    if rawdata.T(i_tend) >= tend
        break;
    end
end

% define trimmed data series
t = rawdata.T(i_tstart : i_tend) - rawdata.T(i_tstart);
x = rawdata.Y(i_tstart : i_tend , :);
u = rawdata.U(i_tstart : i_tend , :);

% save trimmed data files
val_fname = ['trimdataFiles_v2', filesep, data_name, '_trim.mat'];
[unique_fname, change_detect] = auto_rename(val_fname, '0');
save(unique_fname, 't', 'x', 'u');

end
