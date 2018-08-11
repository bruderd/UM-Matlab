function [zmerged, zval, zall] = prep_iddata( data )
%prep_iddata: prepares data to be used with Matlab System Identification
%Toolbox
%   Detailed explanation goes here
%
%   INPUTS:
%       data is a struct containing all of the experimental data for the
%       systm. For more information see functions such as gen_data_fromSim
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
zall = struct;

% calculate the actual total numner of trials
num = round(sqrt(data.valparams.numTrials));
numTrials = num^2;


% initialize merged dataset
zmerged = iddata(data.trial1.y, data.trial1.u, data.valparams.Ts);

% create iddata objects for all trials
for i = 2 : numTrials
   trialID = ['trial', num2str(i)];
   trialName = data.(trialID);
   
   expID = ['z', num2str(i)];
   zall.(expID) = iddata(trialName.y, trialName.u, data.valparams.Ts);
   
   % merge all of the data sets into single multiexperiment object
   zmerged = merge( zmerged, zall.(expID) );
end

% create iddate object for validation data set
zval = iddata(data.validation.x, data.validation.u, data.valparams.Ts);

end

