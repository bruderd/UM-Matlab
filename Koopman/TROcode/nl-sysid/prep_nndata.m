function [ nn_input , nn_output ] = prep_nndata( data )
%prep_nnet: prepares data to be used with Matlab Neural Network Toolbox
%   Detailed explanation goes here
%
%   INPUTS:
%       data is a cell array containing experimental data.
%       It is a cell array where each element is a 
%       struct containing data from a single trial
%
%   OUTPUTS:
%       nn_input is the inputs (u) from all of the trials in data 
%       concatenated together into a single matrix.
%
%       nn_output is the putputs (y) from all of the trials in data 
%       concatenated together into a single matrix.
 
% calculate the actual total numner of trials in data
numTrials = length(data);

% initialize merged dataset
nn_input = data{1}.u;
nn_output = data{1}.y;

% append the inputs/outputs from all trials
for i = 2 : numTrials
   
   nn_input = [ nn_input ; data{i}.u ];
   nn_output = [ nn_output ; data{i}.y ];
   
end

end