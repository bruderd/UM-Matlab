function data = gen_data_fromSim( params )
%genData_fromSim: Generate system data from simulation of known dynamics
%   Detailed explanation goes here

% initialize output struct
data = struct;

%% simulate a bunch of trials with randomized initial conditions

num = round(sqrt(params.numTrials));
ampRange = linspace(params.ampRange(1), params.ampRange(2), num);
freqRange = linspace(params.freqRange(1), params.freqRange(2), num);
trialCount = 1;        % trial counter

x = []; y = [];
for i = 1:num
    for k = 1:num
        
        % generate parameters for input sinusoid
        params.amp = ampRange(i); %5 * (k / num);
        params.freq = freqRange(i); %10 * (i / num);    
    
        % randomize initial conditions
        params.x0 = (params.x0max - params.x0min) .* rand + params.x0min;
            
        % generate data from one simulation
        trialData = run_sim(params);
        
        xk = zeros(length(trialData.t), size(trialData.y,2) + size(trialData.u,2));
        yk = zeros(length(trialData.t), size(trialData.y,2) + size(trialData.u,2));
        for j = 1:length(trialData.t)-1
            xk(j,:) = [ trialData.y(j,:), trialData.u(j,:) ];
            yk(j,:) = [ trialData.y(j+1,:), trialData.u(j,:) ];
        end
        
        % append snapshot pairs from this trial onto set of all pairs
        x = [x; xk];
        y = [y; yk];
        
        % save this trial data to the output struct
        trialID = ['trial', num2str(trialCount)];
        data.(trialID) = trialData;
        
        trialCount = trialCount + 1;    % increment trial counter
    end
end

% define snapshotPairs struct
snapshotPairs = struct;
snapshotPairs.x = x;
snapshotPairs.y = y;

%% Do one simulation to be used for validation (could make this a whole set)

% randomize input and initial contidion
params.amp = (params.ampRange(2) - params.ampRange(1)) .* rand + params.ampRange(1); 
params.freq = (params.freqRange(2) - params.freqRange(1)) .* rand + params.freqRange(1); 
params.x0 = (params.x0max - params.x0min) .* rand + params.x0min;
% params.inputType = 'exponential';

% generate data from one simulation
validation = run_sim(params);


%% Define output

data.snapshotPairs = snapshotPairs;
data.validation = validation;   % trial that can be used for model validation
data.valparams = params;   % saves params used for validation so we can remember

%% save datafile without overwriting previous files with same name
% SaveWithNumber(['dataFiles', filesep, params.systemName, '.mat'], data);
[unique_fname, change_detect] = auto_rename(['dataFiles', filesep, params.systemName, '.mat'], '0');
save(unique_fname, 'data');


end

