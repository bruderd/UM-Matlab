function [data , params] = gen_data_fromExp( params )
%genData_fromExp: Generate system data and snapshot pairs from experimental
%data
%   Detailed explanation goes here

nd = params.nd;  % number of delays included in zeta snapshot pairs

% initialize output struct
data = struct;
alltrials = struct;

%% Read in data from all trials

trialCount = 0;        % trial counter

% Load in sysid data trials
disp('Please select all .mat files corresponding to sysid trials:');
[sysid_file,sysid_path] = uigetfile('MultiSelect','on');
if iscell(sysid_file)
    numTrials = length(sysid_file);
else
    numTrials = 1;
end

alltrials.t = []; alltrials.y = []; alltrials.u = []; alltrials.x = []; alltrials.w = [];
for i = 1 : numTrials
    trialCount = trialCount + 1;    % increment trial counter
    
    % handle single file exception
    if ~iscell(sysid_file)
        file = sysid_file;
    else
        file = sysid_file{i};
    end
    
    % generate data from the ith sysid trial
    trialData = get_data(file, sysid_path, params);
    
    % append this data to the "alltrials" field of data
    alltrials.t = [alltrials.t; trialData.t];     % time vector
    alltrials.y = [alltrials.y; trialData.y];     % state "measurements"
    alltrials.u = [alltrials.u; trialData.u];     % input
    alltrials.x = [alltrials.x; trialData.x];     % actual state
    alltrials.w = [alltrials.w; trialData.w];     % load
    
    % save this trial data to the output struct
    trialID = ['trial', num2str(trialCount)];
    data.(trialID) = trialData;
    
end

%% scale all of the sysid data trials and save snapshot pairs

[alltrials.x, alltrials.u, alltrials.w , params] = scale_data(alltrials.x, alltrials.u, alltrials.w , params);  % scaling factor determined by maxes over all trial data

x = []; y = []; u = []; w = []; zeta_x = []; zeta_y = [];
for i = 1 : numTrials
    
    % rename for convenience
    trialID = ['trial', num2str(i)];
    trialData = data.(trialID);
    
    % scale the data from each of the trials
    trialData.x = trialData.x * diag(params.xScaleFactor);
    trialData.y = trialData.y * diag(params.xScaleFactor);
    trialData.u = trialData.u * diag(params.uScaleFactor);
    trialData.w = trialData.w * diag(params.wScaleFactor);
    
    % pull out snapshot pairs from the scaled data
    n = size(trialData.y , 2);  % length of state vector (also should be specified in params)
    p = size(trialData.u , 2);  % length of input vector (also should be specified in params)
    nw = size(trialData.w , 2);  % length of load vector (also should be specified in params)
    
    xk = zeros(length(trialData.t-1), n);
    yk = zeros(length(trialData.t-1), n);
    wk = zeros(length(trialData.t-1), nw);
    uk = zeros(length(trialData.t-1), p);
    zeta_xk = zeros(length(trialData.t-1), ( size(trialData.y,2) + size(trialData.u,2) ) * nd + size(trialData.y,2) ); % points that include delays
    zeta_yk = zeros(length(trialData.t-1), ( size(trialData.y,2) + size(trialData.u,2) ) * nd + size(trialData.y,2) ); % points that include delays
    for j = nd+1 : length(trialData.t)-1
        xk(j,:) = trialData.y(j,:);
        yk(j,:) = trialData.y(j+1,:);
        uk(j,:) = trialData.u(j,:);
        wk(j,:) = trialData.w(j,:);
        
        % points that include the designated numer of delays
        if nd ~= 0
            zeta_xk(j,:) = [ trialData.y(j,:) ,...
                reshape( flipud( trialData.y(j-nd : j-1, :) )' , [1, n * nd] ) ,...
                reshape( flipud( trialData.u(j-nd : j-1, :) )' , [1, p * nd] ) ];
            zeta_yk(j,:) = [ trialData.y(j+1,:) ,...
                reshape( flipud( trialData.y(j-nd+1 : j, :) )' , [1, n * nd] ) ,...
                reshape( flipud( trialData.u(j-nd+1 : j, :) )' , [1, p * nd] ) ];
        else    % if no delays, zeta_xy is same as xy
            zeta_xk(j,:) = xk(j,:);
            zeta_yk(j,:) = yk(j,:);
        end
    end
    
    % remove initial rows of zeros
    xk = xk(nd+1 : end , :);
    yk = yk(nd+1 : end , :);
    uk = uk(nd+1 : end , :);
    wk = wk(nd+1 : end , :);
    zeta_xk = zeta_xk(nd+1 : end , :);
    zeta_yk = zeta_yk(nd+1 : end , :);
    
    % append snapshot pairs from this trial onto set of all pairs
    x = [x; xk];
    y = [y; yk];
    u = [u; uk];
    w = [w; wk];
    zeta_x = [zeta_x; zeta_xk];
    zeta_y = [zeta_y; zeta_yk];
end

%% define snapshotPairs struct

snapshotPairs = struct;
snapshotPairs.x = x;
snapshotPairs.y = y;
snapshotPairs.u = u;
snapshotPairs.w = w;
snapshotPairs.zeta_x = zeta_x;
snapshotPairs.zeta_y = zeta_y;

%% Read in validation data set(s)

% Load in validation data trials
disp('Please select all .mat files corresponding to validation trials:');
[val_file,val_path] = uigetfile('MultiSelect','on');
if iscell(val_file)
    numVals = length(val_file);
else
    numVals = 1;
end

for j = 1 : numVals
    
    % handle single file exception
    if ~iscell(val_file)
        file = val_file;
    else
        file = val_file{j};
    end
    
    % generate data from the jth validation trial
    valData = get_data(file, val_path, params);
    
    % scale data
    valData.x = valData.x * diag(params.xScaleFactor);
    valData.y = valData.y * diag(params.xScaleFactor);
    valData.u = valData.u * diag(params.uScaleFactor);
    valData.w = valData.w * diag(params.wScaleFactor);
    
    % save this trial data to the output struct
    valID = ['val', num2str(j)];
    data.(valID) = valData;   

end


%% Define output

data.alltrials = alltrials;     % saves data from all trials as a single timeseries
data.snapshotPairs = snapshotPairs;
data.validation = data.val1;   % a trial that can be used for model validation
data.valparams = params;   % saves params used for validation so we can remember

%% save datafile without overwriting previous files with same name
% SaveWithNumber(['dataFiles', filesep, params.systemName, '.mat'], data);
[unique_fname, change_detect] = auto_rename(['dataFiles', filesep, params.system, '.mat'], '0');
save(unique_fname, 'data');


end