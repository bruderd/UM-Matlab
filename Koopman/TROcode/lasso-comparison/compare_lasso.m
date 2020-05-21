% compare_lasso.m
%
% Compare the performance of models across lasso parameters

%% load candidate models
temp = load('..\systems\waves\old-robot_class-nomodels.mat');
class_nomodels = temp.sysid_class;
temp = load('..\systems\waves\old-robot_lasso-candidate-models.mat');
models_proj = temp.candidates;

ksysid_proj = class_nomodels;
ksysid_proj.candidates = models_proj;

%% Run validation on all the models


res_proj = cell( size(ksysid_proj.candidates) );    % store results in a cell array
err_proj = cell( size(ksysid_proj.candidates) );    % store error in a cell array 

if iscell(ksysid.candidates)
    for i = 1 : length(ksysid_proj.candidates)
        [ res_proj{i} , err_proj{i} ] = ksysid_proj.valNplot_model( i );
    end
else
    [ res_proj{1} , err_proj{1} ] = ksysid_proj.valNplot_model;
end