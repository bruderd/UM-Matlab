function [t,x,u] = chop_data(numSlices, rawdata)
%trim_data: trim away the beginning and end of raw data file
%   INPUTS
%       rawdata - struct containing data to be chopped. With fields:
%           .t - time vector
%           .x - state matrix with rows corresponding to state at t(row#)
%           .u - input matrix with rows corresponding to state at t(row#)
%           .y - output matrix with rows corresponding to state at t(row#)
%                (not necessary)
%       numSlices - number of slices the data should be chopped into
%   OUTPUTS
%       t - struct containing time vectors in fields .s1 - .snumSlices
%       x - struct containing state matrix in fields .s1 - .snumSlices
%       u - struct containing input matrix in fields .s1 - .snumSlices

% Check if data argument is provided, it not, promt user to load from file
if ~exist('rawdata', 'var')
    % Prompt user to identify data file
    disp(['Please select .mat file corresponding to raw data.']);
    [data_file,data_path] = uigetfile;
    [data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
    rawdata = load([data_path, data_file]);
end


slice_len = rawdata.t(end) / numSlices; % length of slices, in seconds

last = 0;   % index at which previous slice finished
t = struct;
x = struct;
u = struct;
for i = 1:numSlices
    % find index of first data point past end of slice
    for index = last + 1 : length(rawdata.t)
        if rawdata.t(index) >= i * slice_len
            break;
        end
    end
    
    % construct slice ID
    sliceID = ['s', num2str(i)];
    
    % store ith slice in the ith column of 3rd dimension of t, u, x
    t.(sliceID) = rawdata.t(last + 1 : index) - rawdata.t(last + 1);
    x.(sliceID) = rawdata.x(last + 1 : index , :);
    u.(sliceID) = rawdata.u(last + 1 : index , :);
    
    last = index;
end

 
end