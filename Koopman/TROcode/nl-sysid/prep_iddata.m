function [zsysid_merged, zval_merged, zsysid, zval] = prep_iddata( traindata , valdata )
%prep_iddata: prepares data to be used with Matlab System Identification
%Toolbox
%   Detailed explanation goes here
%
%   INPUTS:
%       traindata/valdata are cell arrays containing all of the experimental
%       data for the systm. It is a cell array where each element is a 
%       struct containing data from a single trial
%
%   OUTPUTS:
%       zmerged is an iddata object containing the merged data from all
%       experements
%
%       zval is an iddata object containing tht data from the validation
%       experiment
%
%       zall is a struct containing the individual data from each
%       experiment. Each experiment can be accessed zall.z#, where # is the
%       experiment's identification number
%
%%
zsysid = struct;
zval = struct;

% calculate the actual total numner of trials
numTrials = length(traindata);
numVals = length(valdata);

%% create iddata objects for sysid trials

% initialize merged dataset
Ts = mean( traindata{1}.t(2:end) - traindata{1}.t(1:end-1) );   % timestep
zsysid.z1 = iddata(traindata{1}.y, traindata{1}.u, Ts);
zsysid_merged = zsysid.z1;

% create iddata objects for all trials
for i = 2 : numTrials
   
   expID = ['z', num2str(i)];
   zsysid.(expID) = iddata(traindata{i}.y, traindata{i}.u, Ts); % assumes all trials have same timestep
   
   % merge all of the data sets into single multiexperiment object
   zsysid_merged = merge( zsysid_merged, zsysid.(expID) );
end

%% create iddata objects for validation trials

% initialize merged dataset
zval.z1 = iddata(valdata{1}.y, valdata{1}.u, Ts);
zval_merged = zval.z1;

% create iddata objects for all trials
for i = 2 : numVals
   
   expID = ['z', num2str(i)];
   zval.(expID) = iddata(valdata{i}.y, valdata{i}.u, Ts);
   
   % merge all of the data sets into single multiexperiment object
   zval_merged = merge( zval_merged, zval.(expID) );
end

% % create iddate object for validation data set
% zval = iddata(data.validation.x, data.validation.u, data.valparams.Ts);

end
